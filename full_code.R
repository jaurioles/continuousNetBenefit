# Define Source Directory
# Change if wanted to run this
setwd("/mnt/bmh01-rds/Sperrin_UKBB_Fairness")
# We also create folders where things within source directory was saved
# As this code was run in a server
if (!dir.exists("Output/study_2_continuous/manuscript/")) {
  dir.create("Output/study_2_continuous/manuscript/", recursive = TRUE)
}
if (!dir.exists("Output/study_2_continuous/manuscript/Plots/")) {
  dir.create("Output/study_2_continuous/manuscript/Plots/", recursive = TRUE)
}



# Set random seed
RANDOM_SEED <- 42
set.seed(RANDOM_SEED)

#------------------#
#### Parameters ####
#------------------#
timepoint <- 3652
number_bootstraps <- 5000 # If reproducing this, you can get some 'close enough' effects with even just 50 bootstraps

# For multiple-treatment cNB
cardiovascular_statins_threshold <- 0.1
cardiovascular_lifestyle_thresholds <- seq(0.01,0.30,0.005)
cardiovascular_lifestyle_weights <- dnorm(cardiovascular_lifestyle_thresholds, mean=0.1,sd=0.02)
cardiovascular_lifestyle_weights[cardiovascular_lifestyle_thresholds>0.1] <- 0
cardiovascular_lifestyle_thresholds <- c(cardiovascular_lifestyle_thresholds, 0.10001)
cardiovascular_lifestyle_weights <- c(cardiovascular_lifestyle_weights, 0)

# For single-treatment varying Statins Thresholds
statinsvar_thresholds <- seq(0.01,0.30,0.005)
# (Weights are multiplied by *(1-statinsvar_thresholds) because of our assumption that d-b is constant)
statinsvar_weights <- dlnorm(statinsvar_thresholds, meanlog = log(0.1), sdlog = 0.3)*(1-statinsvar_thresholds)
plot(statinsvar_thresholds, statinsvar_weights, type = "l", main = "Log-Normal Distribution", xlab = "x", ylab = "Density", xlim = c(0, 0.2))

#-------------------------#
#### Libraries + SetUp ####
#-------------------------#

# Load the required libraries
library(ggplot2)
library(grid)
library(riskCommunicator)
library(survival)
library(mice)
library(tidyverse)
library(riskRegression)
library(WeightedROC) # For weighted AUC
library(pracma)      # For trapezoidal integration

#------------------------------#
#### Functions used in this ####
#------------------------------#

# Creates inverse probability weights for censoring
create_ipwc <- function(df,censoring_formula,timepoint) {
  
  # Fit a logistic regression
  censoring_model <- coxph(censoring_formula,data=df, x = TRUE, y = TRUE)
  
  # What are the weights of the dataframe?
  df$IPWC <- 1/(predictCox(censoring_model,newdata=df,times=timepoint,type="survival")$survival)
  
  # Remove all censored data by time
  df <- df[((df$CENS==0)|(df$TIMECVD>timepoint)),]
  
  return(df)
}

# Calculates net benefit
calcNetBenefit <- function(thresholds,y_true,s_pred,weights=NULL,rescaled=FALSE) {
  
  # If there are no weights given, make them 1
  if (is_empty(weights)) {
    weights <- rep(1,length(y_true))
  }
  
  nb_vector_to_add <- c()
  for (threshold in thresholds) {

    # Threshold prediction for binary
    y_hat <- s_pred>threshold
    
    # Get true positives and false positives vector
    tp_v <- as.numeric((y_true==1)&(y_hat==1))
    fp_v <- as.numeric((y_true==0)&(y_hat==1))
    
    # What is the benefit for each person
    if (rescaled==FALSE) {
      nb_v <-(tp_v - threshold/(1-threshold)*fp_v)*weights
    } else {
      nb_v <-(tp_v/threshold - 1/(1-threshold)*fp_v)*weights
    }
    
    
    # Calculate NB using formula
    nb <- sum(nb_v)/sum(weights)
    
    # Append to vector of all net benefits
    nb_vector_to_add[length(nb_vector_to_add)+1] = nb
    
  }
  return(nb_vector_to_add)
}

# Calculates cardiovascular net benefit
cardNB <- function(y_true,s_pred,weights=NULL) {
  
  # If there are no weights given, make them 1
  if (is_empty(weights)) {
    weights <- rep(1,length(y_true))
  }
  
  # First part of the CV model is to calculate the beenfit of statins
  statins_cnb <- calcNetBenefit(cardiovascular_statins_threshold,y_true,s_pred,weights=weights,rescaled=TRUE)*cardiovascular_statins_threshold
  
  # Second part is to calculate lifestyles average
  lifestyle_nb <- calcNetBenefit(cardiovascular_lifestyle_thresholds,y_true,s_pred,weights=weights,rescaled=TRUE)
  lifestyle_normalisation <- trapz(cardiovascular_lifestyle_thresholds,cardiovascular_lifestyle_weights/cardiovascular_lifestyle_thresholds)
  lifestyle_cnb <- trapz(cardiovascular_lifestyle_thresholds,cardiovascular_lifestyle_weights*lifestyle_nb)/lifestyle_normalisation
  
  # Change normalisation so that they add up to 1
  statins_cnb <- statins_cnb*1
  lifestyle_cnb <- lifestyle_cnb*1
  total_cnb <- statins_cnb+lifestyle_cnb
  
  # Return is list
  cnb <- list()
  cnb$statins <- statins_cnb
  cnb$lifestyle <- lifestyle_cnb
  cnb$total <- total_cnb
  
  return(cnb)
}

# Calculates net benefit for statins with varying threshold
statinsvarNB <- function(y_true,s_pred,weights=NULL) {
  
  # If there are no weights given, make them 1
  if (is_empty(weights)) {
    weights <- rep(1,length(y_true))
  }
  
  # Second part is to calculate statinsvars average
  statinsvar_nb <- calcNetBenefit(statinsvar_thresholds,y_true,s_pred,weights=weights,rescaled=TRUE)
  statinsvar_normalisation <- trapz(statinsvar_thresholds,statinsvar_weights/statinsvar_thresholds)
  statinsvar_cnb <- trapz(statinsvar_thresholds,statinsvar_weights*statinsvar_nb)/statinsvar_normalisation
  
  return(statinsvar_cnb)
}
#---------------------------#
#### Data Pre-Processing ####
#---------------------------#

# Open dataset
cardio_df <- framingham

# Outcome: 10-year cardiovascular death, although we will use survival
cardio_df$y <- as.numeric((framingham$CVD==1)&(framingham$TIMECVD<=3652))

# Is there any previous cardiovascular event?
cardio_df$PREVCVD <- (cardio_df$PREVCHD==1)|(cardio_df$PREVAP==1)|(cardio_df$PREVMI==1)|(cardio_df$PREVSTRK==1)

# Only first observation for each person
cardio_df <- cardio_df[cardio_df$TIME==0,]

# How much missing data
colMeans(is.na(cardio_df)) * 100

# Which covariates to use:
# ESSENTIAL: pspline(AGE, df=3) + SEX + CURSMOKE
# FULLMODEL: pspline(AGE, df=3) + SEX + CURSMOKE + TOTCHOL + SYSBP + BMI + DIABETES + BPMEDS + PREVCVD + PREVHYP 

# BPMEDS PREVCVD PREVHYP
# What is the outcome:
# Surv(TIMECVD, CVD)
formula_cvd_1 <- Y ~ pspline(AGE, df=3) + SEX + CURSMOKE + PREVCVD
formula_cvd_2 <- Y ~ pspline(AGE, df=3) + SEX + CURSMOKE + TOTCHOL + SYSBP + BMI + DIABETES + BPMEDS + PREVCVD + PREVHYP  + GLUCOSE + educ
# This predicts censoring
formula_cens <- Surv(TIMECVD, CENS) ~ pspline(AGE, df=3) + SEX + CURSMOKE + TOTCHOL + SYSBP + BMI + DIABETES + BPMEDS + PREVCVD + PREVHYP  + GLUCOSE + educ

# Only keep relevant data
cardio_df <- cardio_df[c("TIMECVD","CVD", "AGE", "SEX", "CURSMOKE", "TOTCHOL", "SYSBP", 'BMI', 'DIABETES', 'BPMEDS',
                         'PREVCVD', 'PREVHYP','GLUCOSE','educ')]
yy_df <- as.numeric((framingham$CVD==1)&(framingham$TIMECVD<=timepoint))
censored_df <- as.numeric((framingham$CVD==0)&(framingham$TIMECVD<=timepoint))
xx_df <- cardio_df %>% select(-TIMECVD,-CVD)

# Turn relevant variables into factors
cardio_df$SEX <- cardio_df$SEX-1
cardio_df$educ <- as.factor(cardio_df$educ)

# Use MICE to input X
xx_df_mice <- mice(xx_df,maxit = 30, m = 1)
xx_df <- complete(xx_df_mice, action = "stacked")
xx_df <- na.omit(xx_df)

# Put everything back together
cardio_df <- cbind(xx_df,cardio_df %>% select(TIMECVD, CVD))

# Define outcomes and censoring
cardio_df$CENS <- as.integer((cardio_df$CVD==0))
cardio_df$Y <- as.integer((cardio_df$CVD==1)&(cardio_df$TIMECVD<=timepoint))

#----------------------------#
#### Bootstrap Validation ####
#----------------------------#

#### Step 1. Calculate the optimism values.

# For all the bootstraps, the dataset is the og, so weights can already be calculated
bootstrapped_validation <- cardio_df
btpw_val <- create_ipwc(bootstrapped_validation,censoring_formula = formula_cens,timepoint = timepoint)

# Bootstrap - For each bootstrap round, take a sample from original data
# Train model in bootstrapped sample, then test against original data
# Initialise results from bootstraps, including relevant metrics
thresholds <- cardiovascular_lifestyle_thresholds # Thresholds for net benefit

# For Formula 1

# List for correcting optimism
optimism_corrections_log_1 <- list()
optimism_corrections_log_1[["AUROC"]] <- c() # C-statistic
optimism_corrections_log_1[["Slope"]] <- c() # Calibration slope
optimism_corrections_log_1[["CITL"]] <- c()  # Calibration intercept
optimism_corrections_log_1[["OER"]] <- c()   # Observed/Expected ratio
optimism_corrections_log_1[["NB"]] <- list()   # Net Benefit (For thresholds)
optimism_corrections_log_1[["rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
optimism_corrections_log_1[["CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
optimism_corrections_log_1[["CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
optimism_corrections_log_1[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)

# List for correcting mean value and confidence intervals
training_performance_log_1 <- list()
training_performance_log_1[["AUROC"]] <- c() # C-statistic
training_performance_log_1[["Slope"]] <- c() # Calibration slope
training_performance_log_1[["CITL"]] <- c()  # Calibration intercept
training_performance_log_1[["OER"]] <- c()   # Observed/Expected ratio
training_performance_log_1[["NB"]] <- list()   # Net Benefit (For thresholds)
training_performance_log_1[["rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
training_performance_log_1[["TreatAll_NB"]] <- list()   # Net Benefit (For thresholds)
training_performance_log_1[["TreatAll_rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
training_performance_log_1[["TreatNone_NB"]] <- list()   # Net Benefit (For thresholds)
training_performance_log_1[["TreatNone_rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
training_performance_log_1[["CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
training_performance_log_1[["CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_1[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_1[["TreatAll_CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
training_performance_log_1[["TreatAll_CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_1[["TreatAll_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_1[["TreatNone_CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
training_performance_log_1[["TreatNone_CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_1[["TreatNone_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)

# For Formula 2

# List for correcting optimism
optimism_corrections_log_2 <- list()
optimism_corrections_log_2[["AUROC"]] <- c() # C-statistic
optimism_corrections_log_2[["Slope"]] <- c() # Calibration slope
optimism_corrections_log_2[["CITL"]] <- c()  # Calibration intercept
optimism_corrections_log_2[["OER"]] <- c()   # Observed/Expected ratio
optimism_corrections_log_2[["NB"]] <- list()   # Net Benefit (For thresholds)
optimism_corrections_log_2[["rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
optimism_corrections_log_2[["CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
optimism_corrections_log_2[["CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
optimism_corrections_log_2[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)

# List for correcting mean value and confidence intervals
training_performance_log_2 <- list()
training_performance_log_2[["AUROC"]] <- c() # C-statistic
training_performance_log_2[["Slope"]] <- c() # Calibration slope
training_performance_log_2[["CITL"]] <- c()  # Calibration intercept
training_performance_log_2[["OER"]] <- c()   # Observed/Expected ratio
training_performance_log_2[["NB"]] <- list()   # Net Benefit (For thresholds)
training_performance_log_2[["rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
training_performance_log_2[["TreatAll_NB"]] <- list()   # Net Benefit (For thresholds)
training_performance_log_2[["TreatAll_rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
training_performance_log_2[["TreatNone_NB"]] <- list()   # Net Benefit (For thresholds)
training_performance_log_2[["TreatNone_rNB"]] <- list()   # Rescaled Net Benefit (For thresholds)
training_performance_log_2[["CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
training_performance_log_2[["CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_2[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_2[["TreatAll_CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
training_performance_log_2[["TreatAll_CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_2[["TreatAll_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_2[["TreatNone_CNB_Stats"]] <- c() # Continuous Net Benefit (Brought by Statins)
training_performance_log_2[["TreatNone_CNB_Lifes"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_2[["TreatNone_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)


# Start bootstrap for calculating optimism
for (ii in 1:number_bootstraps) {
  
  # Print bootstrap number
  print(ii)
  
  # Sample new dataset
  bootstrapped_indices <- sample(nrow(cardio_df), replace = TRUE)
  bootstrapped_training <- cardio_df[bootstrapped_indices, ]
  
  # Reweight both training and validation samples to account for censoring
  btpw_train <- create_ipwc(bootstrapped_training,censoring_formula = formula_cens,timepoint = timepoint)
  
  # For formula 1
  
  # Fit data
  fit <- glm(formula = formula_cvd_1, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
  
  # Get prediction for training dataset
  s_pred_train <- predict(fit, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element
  # Get prediction for original dataset
  s_pred_val <- predict(fit, newdata = btpw_val, type = "response")
  linear_pred_val <- -1*log((1-s_pred_val)/s_pred_val) # Linear element
  
  # Get metrics and append to list
  
  # AUROC
  auroc_train <- WeightedROC(guess=s_pred_train,label=btpw_train$Y,weight=btpw_train$IPWC) %>% WeightedAUC
  auroc_val <- WeightedROC(guess=s_pred_val,label=btpw_val$Y,weight=btpw_val$IPWC) %>% WeightedAUC
  optimism_corrections_log_1[["AUROC"]][ii] <- auroc_train - auroc_val
  training_performance_log_1[["AUROC"]][[ii]] <- auroc_train
  
  # Calibration - Calculate Slope
  slope_train <- as.numeric(glm(btpw_train$Y~linear_pred_train,family="binomial", weights = btpw_train$IPWC)$coefficients[2])
  slope_val <- as.numeric(glm(btpw_val$Y~linear_pred_val,family="binomial", weights = btpw_val$IPWC)$coefficients[2])
  optimism_corrections_log_1[["Slope"]][ii] <- slope_train - slope_val
  training_performance_log_1[["Slope"]][ii] <- slope_train
  
  # Calibration - Calculate Calibration-On-The-Large
  cotl_train <- as.numeric(glm(btpw_train$Y~offset(linear_pred_train),family="binomial", weights = btpw_train$IPWC)$coefficients)
  cotl_val <- as.numeric(glm(btpw_val$Y~offset(linear_pred_val),family="binomial", weights = btpw_val$IPWC)$coefficients)
  optimism_corrections_log_1[["CITL"]][ii] <- cotl_train - cotl_val
  training_performance_log_1[["CITL"]][ii] <- cotl_train
  
  # Calibration - Observed/Expected ratio
  oer_train <- mean(btpw_train$Y*btpw_train$IPWC)/mean(s_pred_train*btpw_train$IPWC)
  oer_val <- mean(btpw_val$Y*btpw_val$IPWC)/mean(s_pred_val*btpw_val$IPWC)
  optimism_corrections_log_1[["OER"]][ii] <- oer_train - oer_val
  training_performance_log_1[["OER"]][ii] <- oer_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = FALSE)
  nbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val,weights=btpw_val$IPWC,rescaled = FALSE)
  optimism_corrections_log_1[["NB"]][[ii]] <- nbs_train - nbs_val
  training_performance_log_1[["NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = TRUE)
  rnbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val,weights=btpw_val$IPWC,rescaled = TRUE)
  optimism_corrections_log_1[["rNB"]][[ii]] <- rnbs_train - rnbs_val
  training_performance_log_1[["rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  nbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0+1,weights=btpw_val$IPWC,rescaled = FALSE)
  optimism_corrections_log_1[["TreatAll_NB"]][[ii]] <- nbs_train - nbs_val
  training_performance_log_1[["TreatAll_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = TRUE)
  rnbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0+1,weights=btpw_val$IPWC,rescaled = TRUE)
  optimism_corrections_log_1[["TreatAll_rNB"]][[ii]] <- rnbs_train - rnbs_val
  training_performance_log_1[["TreatAll_rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  nbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0+1,weights=btpw_val$IPWC,rescaled = FALSE)
  optimism_corrections_log_1[["TreatNone_NB"]][[ii]] <- nbs_train - nbs_val
  training_performance_log_1[["TreatNone_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC,rescaled = TRUE)
  rnbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0,weights=btpw_val$IPWC,rescaled = TRUE)
  optimism_corrections_log_1[["TreatNone_rNB"]][[ii]] <- rnbs_train - rnbs_val
  training_performance_log_1[["TreatNone_rNB"]][[ii]] <- rnbs_train
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- cardNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  cnb_model_val <- cardNB(btpw_val$Y,s_pred_val,weights=btpw_val$IPWC)
  optimism_corrections_log_1[["CNB_Stats"]][[ii]] <- cnb_model_train$statins - cnb_model_val$statins
  optimism_corrections_log_1[["CNB_Lifes"]][[ii]] <- cnb_model_train$lifestyle - cnb_model_val$lifestyle
  optimism_corrections_log_1[["CNB"]][[ii]] <- cnb_model_train$total - cnb_model_val$total
  training_performance_log_1[["CNB_Stats"]][[ii]] <- cnb_model_train$statins
  training_performance_log_1[["CNB_Lifes"]][[ii]] <- cnb_model_train$lifestyle
  training_performance_log_1[["CNB"]][[ii]] <- cnb_model_train$total
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- cardNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatAll_CNB_Stats"]][[ii]] <- cnb_treatall_train$statins
  training_performance_log_1[["TreatAll_CNB_Lifes"]][[ii]] <- cnb_treatall_train$lifestyle
  training_performance_log_1[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train$total
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- cardNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatNone_CNB_Stats"]][[ii]] <- cnb_treatnone_train$statins
  training_performance_log_1[["TreatNone_CNB_Lifes"]][[ii]] <- cnb_treatnone_train$lifestyle
  training_performance_log_1[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train$total
  
  # For formula 2
  
  # Fit data
  fit <- glm(formula = formula_cvd_2, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
  
  # Get prediction for training dataset
  s_pred_train <- predict(fit, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element
  # Get prediction for original dataset
  s_pred_val <- predict(fit, newdata = btpw_val, type = "response")
  linear_pred_val <- -1*log((1-s_pred_val)/s_pred_val) # Linear element
  
  # Get metrics and append to list
  
  # AUROC
  auroc_train <- WeightedROC(guess=s_pred_train,label=btpw_train$Y,weight=btpw_train$IPWC) %>% WeightedAUC
  auroc_val <- WeightedROC(guess=s_pred_val,label=btpw_val$Y,weight=btpw_val$IPWC) %>% WeightedAUC
  optimism_corrections_log_2[["AUROC"]][ii] <- auroc_train - auroc_val
  training_performance_log_2[["AUROC"]][[ii]] <- auroc_train
  
  # Calibration - Calculate Slope
  slope_train <- as.numeric(glm(btpw_train$Y~linear_pred_train,family="binomial", weights = btpw_train$IPWC)$coefficients[2])
  slope_val <- as.numeric(glm(btpw_val$Y~linear_pred_val,family="binomial", weights = btpw_val$IPWC)$coefficients[2])
  optimism_corrections_log_2[["Slope"]][ii] <- slope_train - slope_val
  training_performance_log_2[["Slope"]][ii] <- slope_train
  
  # Calibration - Calculate Calibration-On-The-Large
  cotl_train <- as.numeric(glm(btpw_train$Y~offset(linear_pred_train),family="binomial", weights = btpw_train$IPWC)$coefficients)
  cotl_val <- as.numeric(glm(btpw_val$Y~offset(linear_pred_val),family="binomial", weights = btpw_val$IPWC)$coefficients)
  optimism_corrections_log_2[["CITL"]][ii] <- cotl_train - cotl_val
  training_performance_log_2[["CITL"]][ii] <- cotl_train
  
  # Calibration - Observed/Expected ratio
  oer_train <- mean(btpw_train$Y*btpw_train$IPWC)/mean(s_pred_train*btpw_train$IPWC)
  oer_val <- mean(btpw_val$Y*btpw_val$IPWC)/mean(s_pred_val*btpw_val$IPWC)
  optimism_corrections_log_2[["OER"]][ii] <- oer_train - oer_val
  training_performance_log_2[["OER"]][ii] <- oer_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = FALSE)
  nbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val,weights=btpw_val$IPWC,rescaled = FALSE)
  optimism_corrections_log_2[["NB"]][[ii]] <- nbs_train - nbs_val
  training_performance_log_2[["NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = TRUE)
  rnbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val,weights=btpw_val$IPWC,rescaled = TRUE)
  optimism_corrections_log_2[["rNB"]][[ii]] <- rnbs_train - rnbs_val
  training_performance_log_2[["rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  nbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0+1,weights=btpw_val$IPWC,rescaled = FALSE)
  optimism_corrections_log_2[["TreatAll_NB"]][[ii]] <- nbs_train - nbs_val
  training_performance_log_2[["TreatAll_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = TRUE)
  rnbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0+1,weights=btpw_val$IPWC,rescaled = TRUE)
  optimism_corrections_log_2[["TreatAll_rNB"]][[ii]] <- rnbs_train - rnbs_val
  training_performance_log_2[["TreatAll_rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  nbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0+1,weights=btpw_val$IPWC,rescaled = FALSE)
  optimism_corrections_log_2[["TreatNone_NB"]][[ii]] <- nbs_train - nbs_val
  training_performance_log_2[["TreatNone_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC,rescaled = TRUE)
  rnbs_val <- calcNetBenefit(thresholds,btpw_val$Y,s_pred_val*0,weights=btpw_val$IPWC,rescaled = TRUE)
  optimism_corrections_log_2[["TreatNone_rNB"]][[ii]] <- rnbs_train - rnbs_val
  training_performance_log_2[["TreatNone_rNB"]][[ii]] <- rnbs_train
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- cardNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  cnb_model_val <- cardNB(btpw_val$Y,s_pred_val,weights=btpw_val$IPWC)
  optimism_corrections_log_2[["CNB_Stats"]][[ii]] <- cnb_model_train$statins - cnb_model_val$statins
  optimism_corrections_log_2[["CNB_Lifes"]][[ii]] <- cnb_model_train$lifestyle - cnb_model_val$lifestyle
  optimism_corrections_log_2[["CNB"]][[ii]] <- cnb_model_train$total - cnb_model_val$total
  training_performance_log_2[["CNB_Stats"]][[ii]] <- cnb_model_train$statins
  training_performance_log_2[["CNB_Lifes"]][[ii]] <- cnb_model_train$lifestyle
  training_performance_log_2[["CNB"]][[ii]] <- cnb_model_train$total
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- cardNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatAll_CNB_Stats"]][[ii]] <- cnb_treatall_train$statins
  training_performance_log_2[["TreatAll_CNB_Lifes"]][[ii]] <- cnb_treatall_train$lifestyle
  training_performance_log_2[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train$total
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- cardNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatNone_CNB_Stats"]][[ii]] <- cnb_treatnone_train$statins
  training_performance_log_2[["TreatNone_CNB_Lifes"]][[ii]] <- cnb_treatnone_train$lifestyle
  training_performance_log_2[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train$total
  
}

# We do a second bootstrap for the training performance CI, for which
# we first train the models in all the data
bootstrapped_training <- cardio_df
# Reweight both training and validation samples to account for censoring
btpw_train <- create_ipwc(bootstrapped_training,censoring_formula = formula_cens,timepoint = timepoint)
fit_1 <- glm(formula = formula_cvd_1, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
fit_2 <- glm(formula = formula_cvd_2, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
# Now calculate performance of same model over bootstraps
for (ii in 1:number_bootstraps) {
  
  # Print bootstrap number
  print(ii)
  
  # Sample new dataset
  bootstrapped_indices <- sample(nrow(cardio_df), replace = TRUE)
  bootstrapped_training <- cardio_df[bootstrapped_indices, ]
  
  # Reweight both training and validation samples to account for censoring
  btpw_train <- create_ipwc(bootstrapped_training,censoring_formula = formula_cens,timepoint = timepoint)
  
  # For formula 1
  
  # Get prediction for training dataset
  s_pred_train <- predict(fit_1, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element
  
  # Get metrics and append to list
  
  # AUROC
  auroc_train <- WeightedROC(guess=s_pred_train,label=btpw_train$Y,weight=btpw_train$IPWC) %>% WeightedAUC
  training_performance_log_1[["AUROC"]][[ii]] <- auroc_train
  
  # Calibration - Calculate Slope
  slope_train <- as.numeric(glm(btpw_train$Y~linear_pred_train,family="binomial", weights = btpw_train$IPWC)$coefficients[2])
  training_performance_log_1[["Slope"]][ii] <- slope_train
  
  # Calibration - Calculate Calibration-On-The-Large
  cotl_train <- as.numeric(glm(btpw_train$Y~offset(linear_pred_train),family="binomial", weights = btpw_train$IPWC)$coefficients)
  training_performance_log_1[["CITL"]][ii] <- cotl_train
  
  # Calibration - Observed/Expected ratio
  oer_train <- mean(btpw_train$Y*btpw_train$IPWC)/mean(s_pred_train*btpw_train$IPWC)
  training_performance_log_1[["OER"]][ii] <- oer_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = FALSE)
  training_performance_log_1[["NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = TRUE)
  training_performance_log_1[["rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  training_performance_log_1[["TreatAll_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = TRUE)
  training_performance_log_1[["TreatAll_rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  training_performance_log_1[["TreatNone_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC,rescaled = TRUE)
  training_performance_log_1[["TreatNone_rNB"]][[ii]] <- rnbs_train
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- cardNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  training_performance_log_1[["CNB_Stats"]][[ii]] <- cnb_model_train$statins
  training_performance_log_1[["CNB_Lifes"]][[ii]] <- cnb_model_train$lifestyle
  training_performance_log_1[["CNB"]][[ii]] <- cnb_model_train$total
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- cardNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train$total
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- cardNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train$total
  
  # For formula 2

  # Get prediction for training dataset
  s_pred_train <- predict(fit_2, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element
  
  # Get metrics and append to list
  
  # AUROC
  auroc_train <- WeightedROC(guess=s_pred_train,label=btpw_train$Y,weight=btpw_train$IPWC) %>% WeightedAUC
  training_performance_log_2[["AUROC"]][[ii]] <- auroc_train
  
  # Calibration - Calculate Slope
  slope_train <- as.numeric(glm(btpw_train$Y~linear_pred_train,family="binomial", weights = btpw_train$IPWC)$coefficients[2])
  training_performance_log_2[["Slope"]][ii] <- slope_train
  
  # Calibration - Calculate Calibration-On-The-Large
  cotl_train <- as.numeric(glm(btpw_train$Y~offset(linear_pred_train),family="binomial", weights = btpw_train$IPWC)$coefficients)
  training_performance_log_2[["CITL"]][ii] <- cotl_train
  
  # Calibration - Observed/Expected ratio
  oer_train <- mean(btpw_train$Y*btpw_train$IPWC)/mean(s_pred_train*btpw_train$IPWC)
  training_performance_log_2[["OER"]][ii] <- oer_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = FALSE)
  training_performance_log_2[["NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = TRUE)
  training_performance_log_2[["rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  training_performance_log_2[["TreatAll_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = TRUE)
  training_performance_log_2[["TreatAll_rNB"]][[ii]] <- rnbs_train
  
  # Net benefit
  nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC,rescaled = FALSE)
  training_performance_log_2[["TreatNone_NB"]][[ii]] <- nbs_train
  
  # Rescaled benefit
  rnbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC,rescaled = TRUE)
  training_performance_log_2[["TreatNone_rNB"]][[ii]] <- rnbs_train
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- cardNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  training_performance_log_2[["CNB_Stats"]][[ii]] <- cnb_model_train$statins
  training_performance_log_2[["CNB_Lifes"]][[ii]] <- cnb_model_train$lifestyle
  training_performance_log_2[["CNB"]][[ii]] <- cnb_model_train$total
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- cardNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatAll_CNB_Stats"]][[ii]] <- cnb_treatall_train$statins
  training_performance_log_2[["TreatAll_CNB_Lifes"]][[ii]] <- cnb_treatall_train$lifestyle
  training_performance_log_2[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train$total
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- cardNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatNone_CNB_Stats"]][[ii]] <- cnb_treatnone_train$statins
  training_performance_log_2[["TreatNone_CNB_Lifes"]][[ii]] <- cnb_treatnone_train$lifestyle
  training_performance_log_2[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train$total
  
}


#### Step 2. Calculate the performance of the model

# For Formula 1
model_performance_1 <- list()

# Thresholds of net benefit
model_performance_1$Threshold <- thresholds

# Remove optimism of values

# AUROC
model_performance_1$AUROC$Val <- (
  mean(unlist(training_performance_log_1$AUROC)) - mean(unlist(optimism_corrections_log_1$AUROC))
)
model_performance_1$AUROC$Low <- (
  quantile(unlist(training_performance_log_1$AUROC),p=0.05) - mean(unlist(optimism_corrections_log_1$AUROC))
)
model_performance_1$AUROC$High <- (
  quantile(unlist(training_performance_log_1$AUROC),p=0.95) - mean(unlist(optimism_corrections_log_1$AUROC))
)

# Slope
# This isn't a great way of calculating this, but not reported in manuscript
model_performance_1$Slope$Val <- (
  mean(unlist(training_performance_log_1$Slope)) - mean(unlist(optimism_corrections_log_1$Slope))
)
model_performance_1$Slope$Low <- (
  quantile(unlist(training_performance_log_1$Slope),p=0.05) - quantile(unlist(optimism_corrections_log_1$Slope),p=0.95)
)
model_performance_1$Slope$High <- (
  quantile(unlist(training_performance_log_1$Slope),p=0.95) - quantile(unlist(optimism_corrections_log_1$Slope),p=0.05)
)

# CITL
# This isn't a great way of calculating this, but not reported in manuscript
model_performance_1$CITL$Val <- (
  mean(unlist(training_performance_log_1$CITL)) - mean(unlist(optimism_corrections_log_1$CITL))
)
model_performance_1$CITL$Low <- (
  quantile(unlist(training_performance_log_1$CITL),p=0.05) - quantile(unlist(optimism_corrections_log_1$CITL),p=0.95)
)
model_performance_1$CITL$High <- (
  quantile(unlist(training_performance_log_1$CITL),p=0.95) - quantile(unlist(optimism_corrections_log_1$CITL),p=0.05)
)

# OER
# This isn't a great way of calculating this, but not reported in manuscript
model_performance_1$OER$Val <- (
  mean(unlist(training_performance_log_1$OER)) - mean(unlist(optimism_corrections_log_1$OER))
)
model_performance_1$OER$Low <- (
  quantile(unlist(training_performance_log_1$OER),p=0.05) - quantile(unlist(optimism_corrections_log_1$OER),p=0.95)
)
model_performance_1$OER$High <- (
  quantile(unlist(training_performance_log_1$OER),p=0.95) - quantile(unlist(optimism_corrections_log_1$OER),p=0.05)
)

# NB
# We don't want the mean here, we calculate again
btpw_train <- create_ipwc(cardio_df,censoring_formula = formula_cens,timepoint = timepoint)
s_pred_train <- predict(fit_1, newdata = btpw_train, type = "response")
nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = FALSE)
model_performance_1$NB$Val <- (
  nbs_train - colMeans(do.call(rbind, optimism_corrections_log_1$NB))
)
model_performance_1$NB$Low <- (
  apply(do.call(rbind, training_performance_log_1$NB), 2, quantile, probs = 0.05,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_1$NB))
)
model_performance_1$NB$High <- (
  apply(do.call(rbind, training_performance_log_1$NB), 2, quantile, probs = 0.95,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_1$NB))
)

# rNB
model_performance_1$rNB$Val <- (
  colMeans(do.call(rbind, training_performance_log_1$rNB)) - colMeans(do.call(rbind, optimism_corrections_log_1$rNB))
)
model_performance_1$rNB$Low <- (
  apply(do.call(rbind, training_performance_log_1$rNB), 2, quantile, probs = 0.05,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_1$rNB))
)
model_performance_1$rNB$High <- (
  apply(do.call(rbind, training_performance_log_1$rNB), 2, quantile, probs = 0.95,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_1$rNB))
)

# TreatAll_NB
# Add
model_performance_1$TreatAll_NB$Val <- (
  colMeans(do.call(rbind, training_performance_log_1$TreatAll_NB))
)
model_performance_1$TreatAll_NB$Low <- (
  apply(do.call(rbind, training_performance_log_1$TreatAll_NB), 2, quantile, probs = 0.05,na.rm=TRUE)
)
model_performance_1$TreatAll_NB$High <- (
  apply(do.call(rbind, training_performance_log_1$TreatAll_NB), 2, quantile, probs = 0.95,na.rm=TRUE)
)

# TreatAll_rNB
model_performance_1$TreatAll_rNB$Val <- (
  colMeans(do.call(rbind, training_performance_log_1$TreatAll_rNB))
)
model_performance_1$TreatAll_rNB$Low <- (
  apply(do.call(rbind, training_performance_log_1$TreatAll_rNB), 2, quantile, probs = 0.05,na.rm=TRUE)
)
model_performance_1$TreatAll_rNB$High <- (
  apply(do.call(rbind, training_performance_log_1$TreatAll_rNB), 2, quantile, probs = 0.95,na.rm=TRUE)
)

# TreatNone_NB
model_performance_1$TreatNone_NB$Val <- (
  thresholds*0
)
model_performance_1$TreatNone_NB$Low <- (
  thresholds*0
)
model_performance_1$TreatNone_NB$High <- (
  thresholds*0
)

# TreatNone_rNB
model_performance_1$TreatNone_rNB$Val <- (
  thresholds*0
)
model_performance_1$TreatNone_rNB$Low <- (
  thresholds*0
)
model_performance_1$TreatNone_rNB$High <- (
  thresholds*0
)

# CNB_Stats
model_performance_1$CNB_Stats$Val <- (
  mean(unlist(training_performance_log_1$CNB_Stats)) - mean(unlist(optimism_corrections_log_1$CNB_Stats))
)
model_performance_1$CNB_Stats$Low <- (
  quantile(unlist(training_performance_log_1$CNB_Stats),p=0.05) - mean(unlist(optimism_corrections_log_1$CNB_Stats))
)
model_performance_1$CNB_Stats$High <- (
  quantile(unlist(training_performance_log_1$CNB_Stats),p=0.95) - mean(unlist(optimism_corrections_log_1$CNB_Stats))
)

# CNB_Lifes
model_performance_1$CNB_Lifes$Val <- (
  mean(unlist(training_performance_log_1$CNB_Lifes)) - mean(unlist(optimism_corrections_log_1$CNB_Lifes))
)
model_performance_1$CNB_Lifes$Low <- (
  quantile(unlist(training_performance_log_1$CNB_Lifes),p=0.05) - mean(unlist(optimism_corrections_log_1$CNB_Lifes))
)
model_performance_1$CNB_Lifes$High <- (
  quantile(unlist(training_performance_log_1$CNB_Lifes),p=0.95) - mean(unlist(optimism_corrections_log_1$CNB_Lifes))
)

# CNB
model_performance_1$CNB$Val <- (
  mean(unlist(training_performance_log_1$CNB)) - mean(unlist(optimism_corrections_log_1$CNB))
)
model_performance_1$CNB$Low <- (
  quantile(unlist(training_performance_log_1$CNB),p=0.05) - mean(unlist(optimism_corrections_log_1$CNB))
)
model_performance_1$CNB$High <- (
  quantile(unlist(training_performance_log_1$CNB),p=0.95) - mean(unlist(optimism_corrections_log_1$CNB))
)

# TreatAll_CNB_Stats
model_performance_1$TreatAll_CNB_Stats$Val <- (
  mean(unlist(training_performance_log_1$TreatAll_CNB_Stats))
)
model_performance_1$TreatAll_CNB_Stats$Low <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB_Stats),p=0.05)
)
model_performance_1$TreatAll_CNB_Stats$High <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB_Stats),p=0.95)
)

# TreatAll_CNB_Lifes
model_performance_1$TreatAll_CNB_Lifes$Val <- (
  mean(unlist(training_performance_log_1$TreatAll_CNB_Lifes))
)
model_performance_1$TreatAll_CNB_Lifes$Low <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB_Lifes),p=0.05)
)
model_performance_1$TreatAll_CNB_Lifes$High <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB_Lifes),p=0.95)
)

# TreatAll_CNB
model_performance_1$TreatAll_CNB$Val <- (
  mean(unlist(training_performance_log_1$TreatAll_CNB))
)
model_performance_1$TreatAll_CNB$Low <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB),p=0.05)
)
model_performance_1$TreatAll_CNB$High <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB),p=0.95)
)

# TreatNone_CNB_Stats
model_performance_1$TreatNone_CNB_Stats$Val <- 0
model_performance_1$TreatNone_CNB_Stats$Low <- 0
model_performance_1$TreatNone_CNB_Stats$High <- 0

# TreatNone_CNB_Lifes
model_performance_1$TreatNone_CNB_Lifes$Val <- 0
model_performance_1$TreatNone_CNB_Lifes$Low <- 0
model_performance_1$TreatNone_CNB_Lifes$High <- 0

# TreatNone_CNB
model_performance_1$TreatNone_CNB$Val <- 0
model_performance_1$TreatNone_CNB$Low <- 0
model_performance_1$TreatNone_CNB$High <- 0





# For Formula 2
model_performance_2 <- list()

# Thresholds of net benefit
model_performance_2$Threshold <- thresholds

# Remove optimism of values

# AUROC
model_performance_2$AUROC$Val <- (
  mean(unlist(training_performance_log_2$AUROC)) - mean(unlist(optimism_corrections_log_2$AUROC))
)
model_performance_2$AUROC$Low <- (
  quantile(unlist(training_performance_log_2$AUROC),p=0.05) - mean(unlist(optimism_corrections_log_2$AUROC))
)
model_performance_2$AUROC$High <- (
  quantile(unlist(training_performance_log_2$AUROC),p=0.95) - mean(unlist(optimism_corrections_log_2$AUROC))
)

# Slope
# This isn't a great way of calculating this, but not reported in manuscript
model_performance_2$Slope$Val <- (
  mean(unlist(training_performance_log_2$Slope)) - mean(unlist(optimism_corrections_log_2$Slope))
)
model_performance_2$Slope$Low <- (
  quantile(unlist(training_performance_log_2$Slope),p=0.05) - quantile(unlist(optimism_corrections_log_2$Slope),p=0.95)
)
model_performance_2$Slope$High <- (
  quantile(unlist(training_performance_log_2$Slope),p=0.95) - quantile(unlist(optimism_corrections_log_2$Slope),p=0.05)
)

# CITL
# This isn't a great way of calculating this, but not reported in manuscript
model_performance_2$CITL$Val <- (
  mean(unlist(training_performance_log_2$CITL)) - mean(unlist(optimism_corrections_log_2$CITL))
)
model_performance_2$CITL$Low <- (
  quantile(unlist(training_performance_log_2$CITL),p=0.05) - quantile(unlist(optimism_corrections_log_2$CITL),p=0.95)
)
model_performance_2$CITL$High <- (
  quantile(unlist(training_performance_log_2$CITL),p=0.95) - quantile(unlist(optimism_corrections_log_2$CITL),p=0.05)
)

# OER
# This isn't a great way of calculating this, but not reported in manuscript
model_performance_2$OER$Val <- (
  mean(unlist(training_performance_log_2$OER)) - mean(unlist(optimism_corrections_log_2$OER))
)
model_performance_2$OER$Low <- (
  quantile(unlist(training_performance_log_2$OER),p=0.05) - quantile(unlist(optimism_corrections_log_2$OER),p=0.95)
)
model_performance_2$OER$High <- (
  quantile(unlist(training_performance_log_2$OER),p=0.95) - quantile(unlist(optimism_corrections_log_2$OER),p=0.05)
)

# NB
# We don't want the mean here, we calculate again
btpw_train <- create_ipwc(cardio_df,censoring_formula = formula_cens,timepoint = timepoint)
s_pred_train <- predict(fit_2, newdata = btpw_train, type = "response")
nbs_train <- calcNetBenefit(thresholds,btpw_train$Y,s_pred_train,weights=btpw_train$IPWC,rescaled = FALSE)
model_performance_2$NB$Val <- (
  nbs_train - colMeans(do.call(rbind, optimism_corrections_log_2$NB))
)
model_performance_2$NB$Low <- (
  apply(do.call(rbind, training_performance_log_2$NB), 2, quantile, probs = 0.05,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_2$NB))
)
model_performance_2$NB$High <- (
  apply(do.call(rbind, training_performance_log_2$NB), 2, quantile, probs = 0.95,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_2$NB))
)

# rNB
model_performance_2$rNB$Val <- (
  colMeans(do.call(rbind, training_performance_log_2$rNB)) - colMeans(do.call(rbind, optimism_corrections_log_2$rNB))
)
model_performance_2$rNB$Low <- (
  apply(do.call(rbind, training_performance_log_2$rNB), 2, quantile, probs = 0.05,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_2$rNB))
)
model_performance_2$rNB$High <- (
  apply(do.call(rbind, training_performance_log_2$rNB), 2, quantile, probs = 0.95,na.rm=TRUE) - colMeans(do.call(rbind, optimism_corrections_log_2$rNB))
)

# TreatAll_NB
model_performance_2$TreatAll_NB$Val <- (
  colMeans(do.call(rbind, training_performance_log_2$TreatAll_NB))
)
model_performance_2$TreatAll_NB$Low <- (
  apply(do.call(rbind, training_performance_log_2$TreatAll_NB), 2, quantile, probs = 0.05,na.rm=TRUE)
)
model_performance_2$TreatAll_NB$High <- (
  apply(do.call(rbind, training_performance_log_2$TreatAll_NB), 2, quantile, probs = 0.95,na.rm=TRUE)
)

# TreatAll_rNB
model_performance_2$TreatAll_rNB$Val <- (
  colMeans(do.call(rbind, training_performance_log_2$TreatAll_rNB))
)
model_performance_2$TreatAll_rNB$Low <- (
  apply(do.call(rbind, training_performance_log_2$TreatAll_rNB), 2, quantile, probs = 0.05,na.rm=TRUE)
)
model_performance_2$TreatAll_rNB$High <- (
  apply(do.call(rbind, training_performance_log_2$TreatAll_rNB), 2, quantile, probs = 0.95,na.rm=TRUE)
)

# TreatNone_NB
model_performance_2$TreatNone_NB$Val <- (
  thresholds*0
)
model_performance_2$TreatNone_NB$Low <- (
  thresholds*0
)
model_performance_2$TreatNone_NB$High <- (
  thresholds*0
)

# TreatNone_rNB
model_performance_2$TreatNone_rNB$Val <- (
  thresholds*0
)
model_performance_2$TreatNone_rNB$Low <- (
  thresholds*0
)
model_performance_2$TreatNone_rNB$High <- (
  thresholds*0
)

# CNB_Stats
model_performance_2$CNB_Stats$Val <- (
  mean(unlist(training_performance_log_2$CNB_Stats)) - mean(unlist(optimism_corrections_log_2$CNB_Stats))
)
model_performance_2$CNB_Stats$Low <- (
  quantile(unlist(training_performance_log_2$CNB_Stats),p=0.05) - mean(unlist(optimism_corrections_log_2$CNB_Stats))
)
model_performance_2$CNB_Stats$High <- (
  quantile(unlist(training_performance_log_2$CNB_Stats),p=0.95) - mean(unlist(optimism_corrections_log_2$CNB_Stats))
)

# CNB_Lifes
model_performance_2$CNB_Lifes$Val <- (
  mean(unlist(training_performance_log_2$CNB_Lifes)) - mean(unlist(optimism_corrections_log_2$CNB_Lifes))
)
model_performance_2$CNB_Lifes$Low <- (
  quantile(unlist(training_performance_log_2$CNB_Lifes),p=0.05) - mean(unlist(optimism_corrections_log_2$CNB_Lifes))
)
model_performance_2$CNB_Lifes$High <- (
  quantile(unlist(training_performance_log_2$CNB_Lifes),p=0.95) - mean(unlist(optimism_corrections_log_2$CNB_Lifes))
)

# CNB
model_performance_2$CNB$Val <- (
  mean(unlist(training_performance_log_2$CNB)) - mean(unlist(optimism_corrections_log_2$CNB))
)
model_performance_2$CNB$Low <- (
  quantile(unlist(training_performance_log_2$CNB),p=0.05) - mean(unlist(optimism_corrections_log_2$CNB))
)
model_performance_2$CNB$High <- (
  quantile(unlist(training_performance_log_2$CNB),p=0.95) - mean(unlist(optimism_corrections_log_2$CNB))
)

# TreatAll_CNB_Stats
model_performance_2$TreatAll_CNB_Stats$Val <- (
  mean(unlist(training_performance_log_2$TreatAll_CNB_Stats))
)
model_performance_2$TreatAll_CNB_Stats$Low <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB_Stats),p=0.05)
)
model_performance_2$TreatAll_CNB_Stats$High <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB_Stats),p=0.95)
)

# TreatAll_CNB_Lifes
model_performance_2$TreatAll_CNB_Lifes$Val <- (
  mean(unlist(training_performance_log_2$TreatAll_CNB_Lifes))
)
model_performance_2$TreatAll_CNB_Lifes$Low <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB_Lifes),p=0.05)
)
model_performance_2$TreatAll_CNB_Lifes$High <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB_Lifes),p=0.95)
)

# TreatAll_CNB
model_performance_2$TreatAll_CNB$Val <- (
  mean(unlist(training_performance_log_2$TreatAll_CNB))
)
model_performance_2$TreatAll_CNB$Low <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB),p=0.05)
)
model_performance_2$TreatAll_CNB$High <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB),p=0.95)
)

# TreatNone_CNB_Stats
model_performance_2$TreatNone_CNB_Stats$Val <- 0
model_performance_2$TreatNone_CNB_Stats$Low <- 0
model_performance_2$TreatNone_CNB_Stats$High <- 0

# TreatNone_CNB_Lifes
model_performance_2$TreatNone_CNB_Lifes$Val <- 0
model_performance_2$TreatNone_CNB_Lifes$Low <- 0
model_performance_2$TreatNone_CNB_Lifes$High <- 0

# TreatNone_CNB
model_performance_2$TreatNone_CNB$Val <- 0
model_performance_2$TreatNone_CNB$Low <- 0
model_performance_2$TreatNone_CNB$High <- 0





# For differences between both formulas
model_performance_diff <- list()

# Thresholds of net benefit
model_performance_diff$Threshold <- thresholds

# Remove optimism of values

# AUROC
model_performance_diff$AUROC$Val <- (
  mean(unlist(training_performance_log_2$AUROC)-unlist(training_performance_log_1$AUROC)) - 
    mean(unlist(optimism_corrections_log_2$AUROC) - unlist(optimism_corrections_log_1$AUROC))
)
model_performance_diff$AUROC$Low <- (
  quantile(unlist(training_performance_log_2$AUROC)-unlist(training_performance_log_1$AUROC),p=0.05) - 
    mean(unlist(optimism_corrections_log_2$AUROC)-unlist(optimism_corrections_log_1$AUROC))
)
model_performance_diff$AUROC$High <- (
  quantile(unlist(training_performance_log_2$AUROC) - unlist(training_performance_log_1$AUROC),p=0.95) - 
    mean(unlist(optimism_corrections_log_2$AUROC) - unlist(optimism_corrections_log_1$AUROC))
)

# CNB_Stats
model_performance_diff$CNB_Stats$Val <- (
  mean(unlist(training_performance_log_2$CNB_Stats) - unlist(training_performance_log_1$CNB_Stats)) - 
    mean(unlist(optimism_corrections_log_2$CNB_Stats) - unlist(optimism_corrections_log_1$CNB_Stats))
)
model_performance_diff$CNB_Stats$Low <- (
  quantile(unlist(training_performance_log_2$CNB_Stats) - unlist(training_performance_log_1$CNB_Stats),p=0.05) - 
    mean(unlist(optimism_corrections_log_2$CNB_Stats) - unlist(optimism_corrections_log_1$CNB_Stats))
)
model_performance_diff$CNB_Stats$High <- (
  quantile(unlist(training_performance_log_2$CNB_Stats) - unlist(training_performance_log_1$CNB_Stats),p=0.95) - 
    mean(unlist(optimism_corrections_log_2$CNB_Stats) - unlist(optimism_corrections_log_1$CNB_Stats))
)

# CNB_Lifes
model_performance_diff$CNB_Lifes$Val <- (
  mean(unlist(training_performance_log_2$CNB_Lifes) - unlist(training_performance_log_1$CNB_Lifes)) - 
    mean(unlist(optimism_corrections_log_2$CNB_Lifes) - unlist(optimism_corrections_log_1$CNB_Lifes))
)
model_performance_diff$CNB_Lifes$Low <- (
  quantile(unlist(training_performance_log_2$CNB_Lifes) - unlist(training_performance_log_1$CNB_Lifes),p=0.05) - 
    mean(unlist(optimism_corrections_log_2$CNB_Lifes) - unlist(optimism_corrections_log_1$CNB_Lifes))
)
model_performance_diff$CNB_Lifes$High <- (
  quantile(unlist(training_performance_log_2$CNB_Lifes) - unlist(training_performance_log_1$CNB_Lifes),p=0.95) - 
    mean(unlist(optimism_corrections_log_2$CNB_Lifes) - unlist(optimism_corrections_log_1$CNB_Lifes))
)

# CNB
model_performance_diff$CNB$Val <- (
  mean(unlist(training_performance_log_2$CNB) - unlist(training_performance_log_1$CNB)) - 
    mean(unlist(optimism_corrections_log_2$CNB) - unlist(optimism_corrections_log_1$CNB))
)
model_performance_diff$CNB$Low <- (
  quantile(unlist(training_performance_log_2$CNB) - unlist(training_performance_log_1$CNB),p=0.05) - 
    mean(unlist(optimism_corrections_log_2$CNB) - unlist(optimism_corrections_log_1$CNB))
)
model_performance_diff$CNB$High <- (
  quantile(unlist(training_performance_log_2$CNB) - unlist(training_performance_log_1$CNB),p=0.95) - 
    mean(unlist(optimism_corrections_log_2$CNB) - unlist(optimism_corrections_log_1$CNB))
)

# Create a model_performance dict for the two formulas and their difference
model_performance <- list()
model_performance$Compact.Model <- model_performance_1
model_performance$Full.Model <- model_performance_2
model_performance$Difference.Models <- model_performance_diff


# Save the results
saveRDS(model_performance,"Output/study_2_continuous/manuscript/cardiovascular_results.rds")

#--------------------------#
#### Validation (Plots) ####
#--------------------------#

# Load the results again
model_performance <- readRDS("Output/study_2_continuous/manuscript/cardiovascular_results.rds")

#------------------------------------#
#### Figure 1a) Different QALYs w ####
#------------------------------------#

# Create data
df <- data.frame(x = seq(0, 15, length.out = 1000))

# Define function
weighting_function <- function(x, a) {
  return(1 / ((1 / x) + (1 / a)))
}

# Add y values for each curve
for (a in seq(0, 1, by = 0.2)) {
  df[paste0("y_a_", a)] <- weighting_function(df$x, a)
}

# Name of label?
lab_name <- "Quality-adjusted\nlife years lost\nwhen treating a\nnegative (d-b)"

# Convert data to long format for ggplot
df_long <- tidyr::pivot_longer(df, cols = starts_with("y_a_"), 
                               names_to = "a",
                               values_to = "y")
df_long$a <- gsub("y_a_", "", df_long$a)

# Plot with ggplot
p <- ggplot(df_long, aes(x = x, y = y, colour = a)) +
  geom_line(size = 1) +
  scale_colour_manual(values = colorRampPalette(c("green", "red"))(length(unique(df_long$a)))) +
  labs(x = 'Quality-adjusted life years added\nwhen treating a positive (a-c)',
       y = 'Weighting',
       colour = lab_name) +
  theme_minimal()

# Save plot
ggsave("Output/study_2_continuous/manuscript/Plots/Fig1a.png", plot = p, width = 5, height = 5, units = "in", dpi = 400)

#------------------------------------#
#### Figure 1b) Different QALYs w ####
#------------------------------------#

# Create data
df <- data.frame(x = seq(0, 1, length.out = 1000))

# Define function
weighting_function <- function(x, a) {
  return(1 / ((1 / x) + (1 / a)))
}

# Add y values for each curve
for (a in seq(0, 15, by = 3)) {
  df[paste0("y_a_", a)] <- weighting_function(df$x, a)
}

# Name of label?
lab_name <- 'Quality-adjusted\nlife years added\nwhen treating a\npositive (a-c)'

# Convert data to long format for ggplot
df_long <- tidyr::pivot_longer(df, cols = starts_with("y_a_"), 
                               names_to = "a",
                               values_to = "y")
df_long$a <- gsub("y_a_", "", df_long$a)
df_long$a <- factor(df_long$a,
                    levels = sort(unique(as.numeric(df_long$a)))
                    )

# Plot with ggplot
p <- ggplot(df_long, aes(x = x, y = y, colour = a)) +
  geom_line(size = 1) +
  scale_colour_manual(values = colorRampPalette(c("green", "red"))(length(unique(df_long$a)))) +
  labs(x = "Quality-adjusted life years lost\nwhen treating a negative (d-b)",
       y = 'Weighting',
       colour = lab_name) +
  theme_minimal()

# Save plot
ggsave("Output/study_2_continuous/manuscript/Plots/Fig1b.png", plot = p, width = 5, height = 5, units = "in", dpi = 400)


#------------------------------#
#### Figure 2a) Net Benefit ####
#------------------------------#

# We first create the dataframe to plot
figure1a_df <- data.frame(
  Optimal.Threshold = model_performance$Compact.Model$Threshold,
  Treat.No.One = 100*model_performance$Compact.Model[["TreatNone_NB"]]$Val,
  Treat.All = 100*model_performance$Compact.Model[["TreatAll_NB"]]$Val,
  Treat.All.LowCI = 100*model_performance$Compact.Model[["TreatAll_NB"]]$Low,
  Treat.All.HighCI = 100*model_performance$Compact.Model[["TreatAll_NB"]]$High,
  Compact.Model = 100*model_performance$Compact.Model[["NB"]]$Val,
  Compact.Model.LowCI = 100*model_performance$Compact.Model[["NB"]]$Low,
  Compact.Model.HighCI = 100*model_performance$Compact.Model[["NB"]]$High,
  Full.Model = 100*model_performance$Full.Model[["NB"]]$Val,
  Full.Model.LowCI = 100*model_performance$Full.Model[["NB"]]$Low,
  Full.Model.HighCI = 100*model_performance$Full.Model[["NB"]]$High
)
figure1a_df <- figure1a_df[-1, ]
figure1a_df <- figure1a_df[-nrow(figure1a_df), ]

# We also need dataframe for lifestyle weight
figure1aweights_df <- data.frame(
  Optimal.Threshold = cardiovascular_lifestyle_thresholds,
  Weights = 0.2*cardiovascular_lifestyle_weights
)

# Plot
p <- ggplot(figure1a_df, aes(x = Optimal.Threshold)) +
  geom_line(aes(y = Treat.No.One, linetype = "Treat.No.One", color = "Treat.No.One")) +
  geom_line(aes(y = Treat.All, linetype = "Treat.All", color = "Treat.All")) +
  geom_line(aes(y = Treat.All.LowCI), linetype = "dotted", color = "black") +
  geom_line(aes(y = Treat.All.HighCI), linetype = "dotted", color = "black") +
  geom_line(aes(y = Compact.Model,linetype = "Compact.Model", color = "Compact.Model")) +
  geom_line(aes(y = Compact.Model.LowCI), linetype = "dotted", color = "blue") +
  geom_line(aes(y = Compact.Model.HighCI), linetype = "dotted", color = "blue") +
  geom_line(aes(y = Full.Model,linetype = "Full.Model", color = "Full.Model")) +
  geom_line(aes(y = Full.Model.LowCI), linetype = "dotted", color = "pink") +
  geom_line(aes(y = Full.Model.HighCI), linetype = "dotted", color = "pink") +
  scale_linetype_manual(values = c("Treat.No.One" = "dotdash",
                                   "Treat.All" = "twodash",
                                   "Compact.Model" = "solid",
                                   "Full.Model" = "solid"),
                        name = "Policy") +
  scale_color_manual(values = c("Treat.No.One" = "grey",
                                   "Treat.All" = "black",
                                "Compact.Model" = "blue",
                                "Full.Model" = "pink"),
                        name = "Policy") +
  labs(x = "Optimal Threshold", y = "Net Benefit (TPs/100 People)") +
  geom_line(data = figure1aweights_df, aes(x=Optimal.Threshold,y=Weights), color ="green") +
  geom_text(data = figure1aweights_df, aes(x = 0.16, y = Weights[1] + 3, label = "Overall Management"), color = "green",size =3) +
  theme_bw() +
  geom_segment(aes(x = 0.1, y = 0, xend = 0.1, yend = 12), arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  geom_text(aes(x = 0.1, y = 12.5, label = "Statins Rx"), vjust = -0.8, color = "red",size =3) +
  coord_cartesian(xlim=c(0,0.2), ylim=c(0,14)) +
  guides(color = guide_legend(title = "Policy"),
         linetype = guide_legend(title = "Policy")) +
  theme(legend.position = c(0.8, 0.8))

# Save plot
ggsave("Output/study_2_continuous/manuscript/Plots/Fig2a.png", plot = p, width = 5, height = 5, units = "in", dpi = 400)

#-----------------------------------------#
#### Figure 2b) Cumulative net benefit ####
#-----------------------------------------#

# Create a dataframe with your data
figure1b_df <- data.frame(
  Model = c("Treat.No.One","Treat.No.One","Treat.No.One",
            "Treat.All","Treat.All","Treat.All",
            "Compact.Model","Compact.Model","Compact.Model",
            "Full.Model","Full.Model","Full.Model"),
  Treatment = c("Statins", "Overall Management","Total Composite",
                "Statins", "Overall Management", "Total Composite",
                "Statins", "Overall Management", "Total Composite",
                "Statins", "Overall Management", "Total Composite"),
  Continuous.Net.Benefit = c(
    100*model_performance$Compact.Model[["TreatNone_CNB_Stats"]]$Val,
    100*model_performance$Compact.Model[["TreatNone_CNB_Lifes"]]$Val,
    100*model_performance$Compact.Model[["TreatNone_CNB"]]$Val,
    100*model_performance$Compact.Model[["TreatAll_CNB_Stats"]]$Val,
    100*model_performance$Compact.Model[["TreatAll_CNB_Lifes"]]$Val,
    100*model_performance$Compact.Model[["TreatAll_CNB"]]$Val,
    100*model_performance$Compact.Model[["CNB_Stats"]]$Val,
    100*model_performance$Compact.Model[["CNB_Lifes"]]$Val,
    100*model_performance$Compact.Model[["CNB"]]$Val,
    100*model_performance$Full.Model[["CNB_Stats"]]$Val,
    100*model_performance$Full.Model[["CNB_Lifes"]]$Val,
    100*model_performance$Full.Model[["CNB"]]$Val
  ),
  Continuous.Net.Benefit.LowCI = c(
    100*model_performance$Compact.Model[["TreatNone_CNB_Stats"]]$Low,
    100*model_performance$Compact.Model[["TreatNone_CNB_Lifes"]]$Low,
    100*model_performance$Compact.Model[["TreatNone_CNB"]]$Low,
    100*model_performance$Compact.Model[["TreatAll_CNB_Stats"]]$Low,
    100*model_performance$Compact.Model[["TreatAll_CNB_Lifes"]]$Low,
    100*model_performance$Compact.Model[["TreatAll_CNB"]]$Low,
    100*model_performance$Compact.Model[["CNB_Stats"]]$Low,
    100*model_performance$Compact.Model[["CNB_Lifes"]]$Low,
    100*model_performance$Compact.Model[["CNB"]]$Low,
    100*model_performance$Full.Model[["CNB_Stats"]]$Low,
    100*model_performance$Full.Model[["CNB_Lifes"]]$Low,
    100*model_performance$Full.Model[["CNB"]]$Low
  ),
  Continuous.Net.Benefit.HighCI = c(
    100*model_performance$Compact.Model[["TreatNone_CNB_Stats"]]$High,
    100*model_performance$Compact.Model[["TreatNone_CNB_Lifes"]]$High,
    100*model_performance$Compact.Model[["TreatNone_CNB"]]$High,
    100*model_performance$Compact.Model[["TreatAll_CNB_Stats"]]$High,
    100*model_performance$Compact.Model[["TreatAll_CNB_Lifes"]]$High,
    100*model_performance$Compact.Model[["TreatAll_CNB"]]$High,
    100*model_performance$Compact.Model[["CNB_Stats"]]$High,
    100*model_performance$Compact.Model[["CNB_Lifes"]]$High,
    100*model_performance$Compact.Model[["CNB"]]$High,
    100*model_performance$Full.Model[["CNB_Stats"]]$High,
    100*model_performance$Full.Model[["CNB_Lifes"]]$High,
    100*model_performance$Full.Model[["CNB"]]$High
  )
)

# Convert Model and Treatment to factor with desired order
figure1b_df$Model <- factor(figure1b_df$Model, levels = c("Treat.No.One", "Treat.All", "Compact.Model", "Full.Model"))
figure1b_df$Treatment <- factor(figure1b_df$Treatment, levels = c("Statins", "Overall Management", "Total Composite"))

# Define custom colours for treatments
treatment_colours <- c("Statins Component" = "red", "Overall Management Component" = "green", "Total Composite" = "purple")
treatment_colours <- c("Treat.No.One" = "grey",
                       "Treat.All" = "black",
                       "Compact.Model" = "blue",
                       "Full.Model" = "pink")

# Remove if total composite
figure1b_df <- figure1b_df[figure1b_df$Treatment!="Total Composite",]

# Plot
p <- ggplot(figure1b_df, aes(x = Treatment, y = Continuous.Net.Benefit, colour = Model, group = Model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Continuous.Net.Benefit.LowCI, ymax = Continuous.Net.Benefit.HighCI), position = position_dodge(width = 0.5), width = 0.2) +
  labs(x = "Intervention/Decision",
       y = "Continuous Net Benefit\n(Statin TPs/100 People)") +
  theme_bw() +
  scale_colour_manual(values = treatment_colours) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Continuous Net Benefit\n(Overall Management TPs/100 People)")) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 8),
        legend.position = "top")

# Save plot
ggsave("Output/study_2_continuous/manuscript/Plots/Fig2b.png", plot = p, width = 5, height = 5, units = "in", dpi = 400)


#### Second Validation (Multiple statins) ####

# In this second validation, we want to calculate the area under the curve
# of statins in a situation where we assume statins side effects are
# constant across the population


#### Step 1. Calculate the optimism values.

# For all the bootstraps, the dataset is the og, so weights can already be calculated
bootstrapped_validation <- cardio_df
btpw_val <- create_ipwc(bootstrapped_validation,censoring_formula = formula_cens,timepoint = timepoint)

# Bootstrap - For each bootstrap round, take a sample from original data
# Train model in bootstrapped sample, then test against original data
# Initialise results from bootstraps, including relevant metrics
thresholds <- cardiovascular_lifestyle_thresholds # Thresholds for net benefit

# For Formula 1

# List for correcting optimism
optimism_corrections_log_1 <- list()
optimism_corrections_log_1[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)

# List for correcting mean value and confidence intervals
training_performance_log_1 <- list()
training_performance_log_1[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_1[["TreatAll_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_1[["TreatNone_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)

# For Formula 2

# List for correcting optimism
optimism_corrections_log_2 <- list()
optimism_corrections_log_2[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)

# List for correcting mean value and confidence intervals
training_performance_log_2 <- list()
training_performance_log_2[["CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_2[["TreatAll_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)
training_performance_log_2[["TreatNone_CNB"]] <- c() # Continuous Net Benefit (Brought by Lifestyle Changes)

for (ii in 1:number_bootstraps) {
  
  # Print bootstrap number
  print(ii)
  
  # Sample new dataset
  bootstrapped_indices <- sample(nrow(cardio_df), replace = TRUE)
  bootstrapped_training <- cardio_df[bootstrapped_indices, ]
  
  # Reweight both training and validation samples to account for censoring
  btpw_train <- create_ipwc(bootstrapped_training,censoring_formula = formula_cens,timepoint = timepoint)
  
  # For formula 1
  
  # Fit data
  fit <- glm(formula = formula_cvd_1, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
  
  # Get prediction for training dataset
  s_pred_train <- predict(fit, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element
  # Get prediction for original dataset
  s_pred_val <- predict(fit, newdata = btpw_val, type = "response")
  linear_pred_val <- -1*log((1-s_pred_val)/s_pred_val) # Linear element
  
  # Get metrics and append to list
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- statinsvarNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  cnb_model_val <- statinsvarNB(btpw_val$Y,s_pred_val,weights=btpw_val$IPWC)
  optimism_corrections_log_1[["CNB"]][[ii]] <- cnb_model_train - cnb_model_val
  training_performance_log_1[["CNB"]][[ii]] <- cnb_model_train
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- statinsvarNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- statinsvarNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train
  
  # For formula 2
  
  # Fit data
  fit <- glm(formula = formula_cvd_2, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
  
  # Get prediction for training dataset
  s_pred_train <- predict(fit, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element
  # Get prediction for original dataset
  s_pred_val <- predict(fit, newdata = btpw_val, type = "response")
  linear_pred_val <- -1*log((1-s_pred_val)/s_pred_val) # Linear element
  
  # Get metrics and append to list
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- statinsvarNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  cnb_model_val <- statinsvarNB(btpw_val$Y,s_pred_val,weights=btpw_val$IPWC)
  optimism_corrections_log_2[["CNB"]][[ii]] <- cnb_model_train - cnb_model_val
  training_performance_log_2[["CNB"]][[ii]] <- cnb_model_train
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- statinsvarNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- statinsvarNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train
  
}

# We do a second bootstrap for the training performance CI, for which
# we first train the models in all the data
bootstrapped_training <- cardio_df
# Reweight both training and validation samples to account for censoring
btpw_train <- create_ipwc(bootstrapped_training,censoring_formula = formula_cens,timepoint = timepoint)
fit_1 <- glm(formula = formula_cvd_1, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
fit_2 <- glm(formula = formula_cvd_2, data = btpw_train, family = binomial(link = logit), weights = btpw_train$IPWC)
# Now calculate performance of same model over bootstraps
for (ii in 1:number_bootstraps) {
  
  # Print bootstrap number
  print(ii)
  
  # Sample new dataset
  bootstrapped_indices <- sample(nrow(cardio_df), replace = TRUE)
  bootstrapped_training <- cardio_df[bootstrapped_indices, ]
  
  # Reweight both training and validation samples to account for censoring
  btpw_train <- create_ipwc(bootstrapped_training,censoring_formula = formula_cens,timepoint = timepoint)
  
  # For formula 1

  # Get prediction for training dataset
  s_pred_train <- predict(fit_1, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element
  
  # Get metrics and append to list
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- statinsvarNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  training_performance_log_1[["CNB"]][[ii]] <- cnb_model_train
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- statinsvarNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- statinsvarNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_1[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train
  
  # For formula 2
  
  # Get prediction for training dataset
  s_pred_train <- predict(fit_2, newdata = btpw_train, type = "response")
  linear_pred_train <- -1*log((1-s_pred_train)/s_pred_train) # Linear element

  # Get metrics and append to list
  
  # Get cardiovascular net benefit for model
  cnb_model_train <- statinsvarNB(btpw_train$Y,s_pred_train,weights=btpw_train$IPWC)
  training_performance_log_2[["CNB"]][[ii]] <- cnb_model_train
  
  # Get cardiovascular net benefit for TreatAll
  cnb_treatall_train <- statinsvarNB(btpw_train$Y,s_pred_train*0+1,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatAll_CNB"]][[ii]] <- cnb_treatall_train
  
  # Get cardiovascular net benefit for TreatNone
  cnb_treatnone_train <- statinsvarNB(btpw_train$Y,s_pred_train*0,weights=btpw_train$IPWC)
  training_performance_log_2[["TreatNone_CNB"]][[ii]] <- cnb_treatnone_train
  
}


#### Step 2. Calculate the performance of the model

# For Formula 1
model_performance_1 <- list()

# Thresholds of net benefit
model_performance_1$Threshold <- thresholds

# Remove optimism of values

# CNB
model_performance_1$CNB$Val <- (
  mean(unlist(training_performance_log_1$CNB)) - mean(unlist(optimism_corrections_log_1$CNB))
)
model_performance_1$CNB$Low <- (
  quantile(unlist(training_performance_log_1$CNB),p=0.05) - mean(unlist(optimism_corrections_log_1$CNB))
)
model_performance_1$CNB$High <- (
  quantile(unlist(training_performance_log_1$CNB),p=0.95) - mean(unlist(optimism_corrections_log_1$CNB))
)

# TreatAll_CNB
model_performance_1$TreatAll_CNB$Val <- (
  mean(unlist(training_performance_log_1$TreatAll_CNB))
)
model_performance_1$TreatAll_CNB$Low <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB),p=0.05)
)
model_performance_1$TreatAll_CNB$High <- (
  quantile(unlist(training_performance_log_1$TreatAll_CNB),p=0.95)
)

# TreatNone_CNB
model_performance_1$TreatNone_CNB$Val <- 0
model_performance_1$TreatNone_CNB$Low <- 0
model_performance_1$TreatNone_CNB$High <- 0





# For Formula 2
model_performance_2 <- list()

# CNB
model_performance_2$CNB$Val <- (
  mean(unlist(training_performance_log_2$CNB)) - mean(unlist(optimism_corrections_log_2$CNB))
)
model_performance_2$CNB$Low <- (
  quantile(unlist(training_performance_log_2$CNB),p=0.05) - mean(unlist(optimism_corrections_log_2$CNB))
)
model_performance_2$CNB$High <- (
  quantile(unlist(training_performance_log_2$CNB),p=0.95) - mean(unlist(optimism_corrections_log_2$CNB))
)


# TreatAll_CNB
model_performance_2$TreatAll_CNB$Val <- (
  mean(unlist(training_performance_log_2$TreatAll_CNB))
)
model_performance_2$TreatAll_CNB$Low <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB),p=0.05)
)
model_performance_2$TreatAll_CNB$High <- (
  quantile(unlist(training_performance_log_2$TreatAll_CNB),p=0.95)
)


# TreatNone_CNB
model_performance_2$TreatNone_CNB$Val <- 0
model_performance_2$TreatNone_CNB$Low <- 0
model_performance_2$TreatNone_CNB$High <- 0





# For differences between both formulas
model_performance_diff <- list()

# Thresholds of net benefit
model_performance_diff$Threshold <- thresholds

# CNB
model_performance_diff$CNB$Val <- (
  mean(unlist(training_performance_log_2$CNB) - unlist(training_performance_log_1$CNB)) - 
    mean(unlist(optimism_corrections_log_2$CNB) - unlist(optimism_corrections_log_1$CNB))
)
model_performance_diff$CNB$Low <- (
  quantile(unlist(training_performance_log_2$CNB) - unlist(training_performance_log_1$CNB),p=0.05) - 
    mean(unlist(optimism_corrections_log_2$CNB) - unlist(optimism_corrections_log_1$CNB))
)
model_performance_diff$CNB$High <- (
  quantile(unlist(training_performance_log_2$CNB) - unlist(training_performance_log_1$CNB),p=0.95) - 
    mean(unlist(optimism_corrections_log_2$CNB) - unlist(optimism_corrections_log_1$CNB))
)

# Create a model_performance dict for the two formulas and their difference
model_performance <- list()
model_performance$Compact.Model <- model_performance_1
model_performance$Full.Model <- model_performance_2
model_performance$Difference.Models <- model_performance_diff

# Save the results
saveRDS(model_performance,"Output/study_2_continuous/manuscript/varstatins_results.rds")


#### Second validation (Plots) ####

# Load the results again
model_performance <- readRDS("Output/study_2_continuous/manuscript/cardiovascular_results.rds")
statvars_performance <- readRDS("Output/study_2_continuous/manuscript/varstatins_results.rds")

#### Figure 3a) Net Benefit ####
#------------------------------#

# We first create the dataframe to plot
figure3a_df <- data.frame(
  Optimal.Threshold = model_performance$Compact.Model$Threshold,
  Treat.No.One = 100*model_performance$Compact.Model[["TreatNone_NB"]]$Val,
  Treat.All = 100*model_performance$Compact.Model[["TreatAll_NB"]]$Val,
  Treat.All.LowCI = 100*model_performance$Compact.Model[["TreatAll_NB"]]$Low,
  Treat.All.HighCI = 100*model_performance$Compact.Model[["TreatAll_NB"]]$High,
  Compact.Model = 100*model_performance$Compact.Model[["NB"]]$Val,
  Compact.Model.LowCI = 100*model_performance$Compact.Model[["NB"]]$Low,
  Compact.Model.HighCI = 100*model_performance$Compact.Model[["NB"]]$High,
  Full.Model = 100*model_performance$Full.Model[["NB"]]$Val,
  Full.Model.LowCI = 100*model_performance$Full.Model[["NB"]]$Low,
  Full.Model.HighCI = 100*model_performance$Full.Model[["NB"]]$High
)
figure3a_df <- figure3a_df[-1, ]
figure3a_df <- figure3a_df[-nrow(figure3a_df), ]

# We also need dataframe for lifestyle weight
figure1aweights_df <- data.frame(
  Optimal.Threshold = statinsvar_thresholds,
  Weights = 0.2*statinsvar_weights
)

# Plot
p <- ggplot(figure3a_df, aes(x = Optimal.Threshold)) +
  geom_line(aes(y = Treat.No.One, linetype = "Treat.No.One", color = "Treat.No.One")) +
  geom_line(aes(y = Treat.All, linetype = "Treat.All", color = "Treat.All")) +
  geom_line(aes(y = Treat.All.LowCI), linetype = "dotted", color = "black") +
  geom_line(aes(y = Treat.All.HighCI), linetype = "dotted", color = "black") +
  geom_line(aes(y = Compact.Model,linetype = "Compact.Model", color = "Compact.Model")) +
  geom_line(aes(y = Compact.Model.LowCI), linetype = "dotted", color = "blue") +
  geom_line(aes(y = Compact.Model.HighCI), linetype = "dotted", color = "blue") +
  geom_line(aes(y = Full.Model,linetype = "Full.Model", color = "Full.Model")) +
  geom_line(aes(y = Full.Model.LowCI), linetype = "dotted", color = "pink") +
  geom_line(aes(y = Full.Model.HighCI), linetype = "dotted", color = "pink") +
  scale_linetype_manual(values = c("Treat.No.One" = "dotdash",
                                   "Treat.All" = "twodash",
                                   "Compact.Model" = "solid",
                                   "Full.Model" = "solid"),
                        name = "Policy") +
  scale_color_manual(values = c("Treat.No.One" = "grey",
                                "Treat.All" = "black",
                                "Compact.Model" = "blue",
                                "Full.Model" = "pink"),
                     name = "Policy") +
  labs(x = "Optimal Threshold", y = "Net Benefit (TPs/100 People)") +
  geom_line(data = figure1aweights_df, aes(x=Optimal.Threshold,y=Weights), color ="red") +
  geom_text(data = figure1aweights_df, aes(x = 0.05, y = 3, label = "Statins\nDistribution"), color = "red",size =3) +
  theme_bw() +
  coord_cartesian(xlim=c(0,0.2), ylim=c(0,14)) +
  guides(color = guide_legend(title = "Policy"),
         linetype = guide_legend(title = "Policy")) +
  theme(legend.position = c(0.8, 0.8))

# Save plot
ggsave("Output/study_2_continuous/manuscript/Plots/Fig3a.png", plot = p, width = 5, height = 5, units = "in", dpi = 400)

#-----------------------------------------#
#### Figure 3b) Cumulative vs Pointwise Net Benefit ####
#-----------------------------------------#

# Create a dataframe with your data
figure1b_df <- data.frame(
  Model = c("Treat.No.One","Treat.No.One","Treat.No.One",
            "Treat.All","Treat.All","Treat.All",
            "Compact.Model","Compact.Model","Compact.Model",
            "Full.Model","Full.Model","Full.Model"),
  Treatment = c("Statins\n(10%)", "Statins\n(Weighted)","Total Composite",
                "Statins\n(10%)", "Statins\n(Weighted)", "Total Composite",
                "Statins\n(10%)", "Statins\n(Weighted)", "Total Composite",
                "Statins\n(10%)", "Statins\n(Weighted)", "Total Composite"),
  Continuous.Net.Benefit = c(
    100*model_performance$Compact.Model[["TreatNone_CNB_Stats"]]$Val,
    100*statvars_performance$Compact.Model[["TreatNone_CNB"]]$Val,
    100*model_performance$Compact.Model[["TreatNone_CNB"]]$Val,
    100*model_performance$Compact.Model[["TreatAll_CNB_Stats"]]$Val,
    100*statvars_performance$Compact.Model[["TreatAll_CNB"]]$Val,
    100*model_performance$Compact.Model[["TreatAll_CNB"]]$Val,
    100*model_performance$Compact.Model[["CNB_Stats"]]$Val,
    100*statvars_performance$Compact.Model[["CNB"]]$Val,
    100*model_performance$Compact.Model[["CNB"]]$Val,
    100*model_performance$Full.Model[["CNB_Stats"]]$Val,
    100*statvars_performance$Full.Model[["CNB"]]$Val,
    100*model_performance$Full.Model[["CNB"]]$Val
  ),
  Continuous.Net.Benefit.LowCI = c(
    100*model_performance$Compact.Model[["TreatNone_CNB_Stats"]]$Low,
    100*statvars_performance$Compact.Model[["TreatNone_CNB"]]$Low,
    100*model_performance$Compact.Model[["TreatNone_CNB"]]$Low,
    100*model_performance$Compact.Model[["TreatAll_CNB_Stats"]]$Low,
    100*statvars_performance$Compact.Model[["TreatAll_CNB"]]$Low,
    100*model_performance$Compact.Model[["TreatAll_CNB"]]$Low,
    100*model_performance$Compact.Model[["CNB_Stats"]]$Low,
    100*statvars_performance$Compact.Model[["CNB"]]$Low,
    100*model_performance$Compact.Model[["CNB"]]$Low,
    100*model_performance$Full.Model[["CNB_Stats"]]$Low,
    100*statvars_performance$Full.Model[["CNB"]]$Low,
    100*model_performance$Full.Model[["CNB"]]$Low
  ),
  Continuous.Net.Benefit.HighCI = c(
    100*model_performance$Compact.Model[["TreatNone_CNB_Stats"]]$High,
    100*statvars_performance$Compact.Model[["TreatNone_CNB"]]$High,
    100*model_performance$Compact.Model[["TreatNone_CNB"]]$High,
    100*model_performance$Compact.Model[["TreatAll_CNB_Stats"]]$High,
    100*statvars_performance$Compact.Model[["TreatAll_CNB"]]$High,
    100*model_performance$Compact.Model[["TreatAll_CNB"]]$High,
    100*model_performance$Compact.Model[["CNB_Stats"]]$High,
    100*statvars_performance$Compact.Model[["CNB"]]$High,
    100*model_performance$Compact.Model[["CNB"]]$High,
    100*model_performance$Full.Model[["CNB_Stats"]]$High,
    100*statvars_performance$Full.Model[["CNB"]]$High,
    100*model_performance$Full.Model[["CNB"]]$High
  )
)

# Convert Model and Treatment to factor with desired order
figure1b_df$Model <- factor(figure1b_df$Model, levels = c("Treat.No.One", "Treat.All", "Compact.Model", "Full.Model"))
figure1b_df$Treatment <- factor(figure1b_df$Treatment, levels = c("Statins\n(10%)", "Statins\n(Weighted)", "Total Composite"))

# Define custom colours for treatments
treatment_colours <- c("Statins Component" = "red", "Overall Management Component" = "green", "Total Composite" = "purple")
treatment_colours <- c("Treat.No.One" = "grey",
                       "Treat.All" = "black",
                       "Compact.Model" = "blue",
                       "Full.Model" = "pink")

# Remove if total composite
figure1b_df <- figure1b_df[figure1b_df$Treatment!="Total Composite",]

# Plot
p <- ggplot(figure1b_df, aes(x = Treatment, y = Continuous.Net.Benefit, colour = Model, group = Model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Continuous.Net.Benefit.LowCI, ymax = Continuous.Net.Benefit.HighCI), position = position_dodge(width = 0.5), width = 0.2) +
  labs(x = "Optimal Threshold Assumption",
       y = "Continuous Net Benefit\n(TPs/100 People)") +
  theme_bw() +
  scale_colour_manual(values = treatment_colours) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Continuous Net Benefit\n(Average TPs/100 People)")) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 8),
        legend.position = "top")

# Save plot
ggsave("Output/study_2_continuous/manuscript/Plots/Fig3b.png", plot = p, width = 5, height = 5, units = "in", dpi = 400)

