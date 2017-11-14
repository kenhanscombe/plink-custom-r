
Rplink <- function(PHENO, GENO, CLUSTER, COVAR) {
  
  library(tidyverse)
  library(broom)
  
  pseudo_rsq <- function(model){
    dev <- model$deviance
    null_dev <- model$null.deviance
    model_n <- length(model$fitted.values)
    r2_cox_snell <- 1 - exp(-(null_dev - dev) / model_n)
    r2_nagelkerke <- r2_cox_snell / (1 - (exp(-(null_dev / model_n))))
    r2_nagelkerke
  }
  
  func <- function(snp) { 
    m <- glm(PHENO == 2 ~ COVAR + snp, family = "binomial")
    rsq <- pseudo_rsq(m)
    glance_m  <- glance(m) %>% unlist(.[1, ])
    tidy_m <- tidy(m) %>% select(-term) %>% tail(n = 1) %>% unlist()
    summary_m <- c(tidy_m, glance_m, rsq)
    c(length(summary_m), summary_m)
  }
  
  apply(GENO, 2, func)
}
