library(boot)
############## somewhat general TMLE framework takes initial data, an estimation function and an update
############## function can be used for arbitrary TMLEs

# truncation function for Q so logistic regression doesn't break on Y close to 0 or 1
#' @export
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
    pmin(pmax(x, lower), upper)
}

# function to estimate logistic parametric submodel and get updated estimate logistic
# fluctuation
#' @export
logit_fluctuate <- function(tmledata, flucmod, truncate = 0) {
    suppressWarnings({
        fluc <- glm(flucmod, data = tmledata, family = "binomial")
    })
    list(eps = coef(fluc))
}


#' @title gentmle
#' @description General TMLE function that takes care of the bookkeeping of estimation and update steps.
#' 
#' @param estimate_fun Function for estimation step
#' @param update_fun, Function for update step
#' @param max_iter, Maximum number of iteration steps
#' @param ..., Extra arguments that can be passed to update_fun and estimate_fun
#'
#' @export
#' @example /inst/examples/ey1_example.R
#' @example /inst/examples/eysi_example.R
gentmle1 <- function(initdata, estimate_fun, update_fun, max_iter = 100, ...) {
    converge <- F
    n=nrow(initdata)
    # cat(sprintf('bw: %f\n',bw))
    eststep <- estimate_fun(initdata)
    
    initests <- eststep$ests
    order <- 1/nrow(initdata)
    e=0
    for (j in seq_len(max_iter)) {
      
      if ((sum(as.numeric(eststep$tmledata$HA==0))==n)|
        (sum(as.numeric(eststep$tmledata$CYA==0))==n))
      {
        ED <- sapply(eststep$Dstar, mean)
        break}
        updatestep <- update_fun(eststep$tmledata, ...)
        e = updatestep$coefs[[1]]
        eststep <- estimate_fun(updatestep$tmledata, ...)
        
        ED <- sapply(eststep$Dstar, mean)
        # cat(sprintf('ED_psi=%e ED_sigma=%e psi=%f sigma2=%f\n coef_h=%f coef_Cy=%f
        # coef_Cg=%f\n',ED[1],ED[2],eststep$ests[1],sigma2=eststep$ests[2],updatestep$coefs[1],updatestep$coefs[2],updatestep$coefs[3]))
        
        if (all(abs(ED) < order)) {
            converge <- T
            break
        }
        
        
        
    }
    
    ED2 <- sapply(eststep$Dstar, function(x) mean(x^2))
    ED3 <- sapply(eststep$Dstar, function(x) mean(x^3))
    result <- list(initdata = initdata, tmledata = eststep$tmledata, initests = initests, tmleests = eststep$ests, 
        steps = j, Dstar = eststep$Dstar, ED = ED, ED2 = ED2, ED3 = ED3,e=e,
        Q1star=eststep$tmledata$Q1k)
    
    # class(result) <- "gentmle"
    
    return(result)
}


#' @export
ci_gentmle1 <- function(gentmle_obj, level = 0.95) {
    
    n <- nrow(gentmle_obj$initdata)
    n_ests <- length(gentmle_obj$tmleests)
    ldply(seq_len(n_ests), function(i) {
        
        est <- gentmle_obj$tmleests[i]
        sd <- sqrt(gentmle_obj$ED2[i])/sqrt(n)
        z <- (1 + level)/2
        lower <- est - qnorm(z) * sd
        upper <- est + qnorm(z) * sd
        e = gentmle_obj$e
        data.frame(parameter = names(est), est = est, sd = sd, lower = lower, upper = upper,e=e)
    })
}

#' @export
print.gentmle1 <- function(gentmle_obj) {
    cat(sprintf("TMLE ran for %d step(s)\n", gentmle_obj$steps))
    EDtext <- sprintf("E[%s]=%1.2e", names(gentmle_obj$ED), gentmle_obj$ED)
    cat(sprintf("The mean of the IC is %s\n", paste(EDtext, collapse = ", ")))
    
    cat("\n\n")
    print(ci_gentmle(gentmle_obj))
    
    cat("\n")
    
} 
