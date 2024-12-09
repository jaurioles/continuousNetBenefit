# Import library
library(pracma) # For trapezoidal integration
library(riskCommunicator) # To import example

#### Calculates the net benefit of a model for: ####
#   - true outcomes of individuals `y_true`
#   - predictions of the model `s_pred`
#   - a set of optimal thresholds `thresholds`
#   - whether the output of the function is the regular
#     net benefit (`rescaled`=FALSE), or the rescaled version
#     (`rescaled`=TRUE)
calcNB <- function(y_true, s_pred, thresholds, rescaled=FALSE) {
  
  # Calculate the net benefit for each threshold on loop
  nb_vector_to_add <- c()
  for (threshold in thresholds) {
    
    # Binarise predictions intro treat/no-treat
    y_hat <- s_pred>threshold
    
    # Get true positives and false positives vector
    tp_v <- as.numeric((y_true==1)&(y_hat==1))
    fp_v <- as.numeric((y_true==0)&(y_hat==1))
    
    # Calculate the net benefit in each individual,
    if (rescaled==FALSE) {
      nb_v <-(tp_v - threshold/(1-threshold)*fp_v)
    } else {
      nb_v <- (tp_v/threshold - fp_v/(1-threshold))
    }
    
    # Average over sample
    nb <- sum(nb_v)/length(nb_v)
    
    # Append to vector of all net benefits
    nb_vector_to_add[length(nb_vector_to_add)+1] = nb
    
  }
  
  # Return list of nb
  return(nb_vector_to_add)
}




#### Calculates continuous net benefit for: ####
#   - true outcomes of individuals `y_true`
#   - predictions of the model `s_pred`
#   - a range of thresholds in which to evaluate the integral `thresholds`
#   - the importance function across these thresholds `importance`
calcContinuousNB <- function(y_true, s_pred, thresholds, importance) {
  
  # Second part is to calculate lifestyles average
  nb <- calcNB(y_true,s_pred,thresholds,rescaled=TRUE)

  # Integrates over range to calculate the continuous net benefit  
  cnb <- trapz(thresholds,nb*importance)
  
  # Normalises the cnb to make sure that the final unit is true positives
  normalisation_scaling <- trapz(thresholds,importance/thresholds)
  cnb <- cnb/normalisation_scaling
  
  
  return(cnb)
}




#### Example of use ####

# Import example dataset
df <- framingham

# Outcome: 10-year cardiovascular death, although we will use survival
df$y <- as.numeric((framingham$CVD==1)&(framingham$TIMECVD<=3652))

# Only first observation for each person
df <- df[df$TIME==0,]

# Train simple logistic regression
fit <- glm(formula = y ~ AGE + SEX + CURSMOKE, data = df, family = binomial(link = logit))

# Create predictions of model
preds <- predict(fit, newdata = df, type = "response")

# To calculate the continuous net benefit, we choose as an example
# weighting function a gaussian centered at 10%
example_thresholds <- seq(0.001,0.3,0.001) # Range of thresholds where importance is non-negligible
example_importance <- dnorm(example_thresholds, mean=0.1,sd=0.02) # Mean 10%, SD 2%

# Calculate the continuous net benefit
cnb <- calcContinuousNB(df$y, preds, example_thresholds, example_importance)

# Print
print(sprintf("The continuous net benefit of the model is %.1f TPs per 100 people.",cnb*100))