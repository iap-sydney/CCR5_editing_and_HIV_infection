
# -------------------------------------------------------------------------
#' Steffen Docken
#' 20-3-23
#' corrected or uncorrected AIC calculation
#' 
# -------------------------------------------------------------------------

# lnL: maximum log-likelihood 
# k: number of parameters
# n: number of measurements
# Corrected: 0 (uncorrected AIC), 1 (corrected AIC)

AIC_generic <- function(lnL, k, n, Corrected){ 
  
  AIC_uncorrect = 2*k - 2*lnL;
  
  if (Corrected == 0){
    AIC_val = AIC_uncorrect;
  }else{
    AIC_val = AIC_uncorrect + (2*k^2 +2*k)/(n - k - 1);
  }
  
  return(AIC_val)
}