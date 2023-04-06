##' Function crossvalid
##'
##' @title Cross-validation schemes
##' 
##' @description
##' This functions performs cross-validations for obtaining the predictive ability
##' of a genomic prediction model. The function can perform CV1, CV2, CV0 and CV00.
##' The models are fitted using [BGLR::BGLR()].
##' 
##' @details
##' Cross-validations are useful for assessing how good the model can be in predicting
##' the performance of genotypes in distinct scenarios. Each cross-validation 
##' scheme (CV1, CV2, CV0 and CV00) mimics a different situation. At CV1, the goal
##' is to predict the performance of untested genotypes. At CV2, genotypes are 
##' evaluated at some environments, but not all of them, so the goal is to predict
##' the performance of genotypes in environment where it was not tested. CV0 assess 
##' the model's ability to predict the performance of tested genotypes in untested 
##' environments. Finally, CV00 evaluates how good is the model's predictions for 
##' untested genotypes performance in untested environments.
##' 
##' @param data A data frame containing the observed data
##' @param y The name of the column that contains the trait values
##' @param gen The name of the column that contains the genotypes information
##' @param env The name of the column that contains the environments information
##' @param ETA A two-level list used to specify the regression function (or linear 
##' predictor). See [BGLR::BGLR()] manual for more details.
##' @param niter An integer indicating the number of iteration
##' @param burnin An integer indicating the number of burn-in iterations
##' @param seed An integer indicating the seed to be used. Mandatory, for the sake
##' of reproducibility
##' @param nfolds An integer indicating the number of folds for the cross-validation
##' @param nrept An integer indicating the number of times the cross-validation 
##' should repeat
##' @param cv 1 for CV1 (default), 2 for CV2, 0 for CV0, and 00 for CV00
##' @param results If NULL, the results will be provided per repetition. 
##' If `results = 'mean'`, the function will return the mean across repetitions.
##' 
##' @return The function will return a list with:
##' \itemize{
##' \item \code{yhat} : A data frame containing the observed and predicted values
##' for each genotype, in each environment.
##' \item \code{corr} : The correlation between observed and predicted values.
##' \item \code{corr_env} : The correlation between observed and predicted values
##'  per environment
##' \item \code{wapa} : The weighted average predictive ability of the model
##' \item \code{mspe} : The mean square prediction error of the model
##' \item \code{slope} : The value of the slope of a linear regression between the 
##' observed and predicted values
##' } 
##' 
##' @import BGLR lubridate
##' 
##' @author Saulo F. S. Chaves (saulo.chaves at ufv.br)
##'

crossvalid = function(data, y, gen, env, ETA, niter, burnin, seed, 
                      nfolds, nrept, cv = 1, results = NULL,...){
  
  data = data
  data$gen = data[,gen]
  data$env = data[,env]
  
  pb <- txtProgressBar(min = 0,
                       max = nrept,
                       style = 3,
                       width = nrept, 
                       char = "=") 
  
  init <- numeric(nrept)
  end <- numeric(nrept)
  
  if (cv == 1){
    
    set.seed(seed)
    sets = split(
      rep(
        1:nfolds, length(unique(data$gen)) * nrept
      )[order(runif(length(unique(data$gen)) * nrept))],
      f = 1:nrept
    )
    reps = list()
    for (j in 1:nrept) {
      init[j] = Sys.time()
      means = data[,c('gen','env',y)]
      means = merge(means, data.frame(
        gen = unique(data$gen),
        set = sets[[j]]
      ), by = gen)
      means$Yhat = NA
      for (i in 1:nfolds) {
        
        yNA = means[,y]
        yNA[means$set == i] = NA
        
        crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = niter, burnIn = burnin, 
                              verbose = F, ...)
        unlink("*.dat")
        
        means$Yhat[means$set == i] = crossval$yHat[means$set == i]
      }
      reps[[j]] = means
      rm(means, crossval, yNA)
      end[j] = Sys.time()
      setTxtProgressBar(pb, j)
      time = round(lubridate::seconds_to_period(sum(end-init)),0)
      est = nrept * (mean(end[end != 0] - init[init != 0])) - time
      remaining = round(lubridate::seconds_to_period(est), 0)
      
      cat(paste(" // Execution time:", time,
                " // Time you have to make a coffee", remaining), "")
    }
    close(pb)
    corr = lapply(reps, function(x) cor(x[, y], x[,'Yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c(y, 'Yhat')], x[,'env'], 
                   function(z) cor(z[,y], z[,'Yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(data[, 'gen'], data[,'env'], length)-2)))/
        (sum(1/(tapply(data[, 'gen'], data[, 'env'], length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'Yhat'] - x[, y])^2))
    slope = lapply(reps, function(x) unname(lm(x[,y] ~ x[,'Yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps_M1_CV1,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }else if (cv == 2){
    reps = list()
    for (j in 1:nrept) {
      init[j] = Sys.time()
      set.seed(seed+j)
      means = data[,c('gen','env',y)]
      means$set = NA
      for (i in unique(data$gen)) {
        means[means$gen == i,'set'] = sample(
          1:nfolds, 
          size = dim(means[means$gen == i,])[1],
          replace = dim(means[means$gen == i,])[1] > nfolds
        )
      }
      means$Yhat = NA
      for (i in 1:nfolds) {
        
        yNA = means[,y]
        yNA[means$set == i] = NA
        
        crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = niter, burnIn = burnin, 
                              verbose = F, ...)
        unlink("*.dat")
        
        means$Yhat[means$set == i] = crossval$yHat[means$set == i]
      }
      reps[[j]] = means
      rm(means, crossval, yNA)
      end[j] = Sys.time()
      setTxtProgressBar(pb, j)
      time = round(lubridate::seconds_to_period(sum(end-init)),0)
      est = nrept * (mean(end[end != 0] - init[init != 0])) - time
      remaining = round(lubridate::seconds_to_period(est), 0)
      cat(paste(" // Execution time:", time,
                " // Time you have to make a coffee", remaining), "")
    }
    corr = lapply(reps, function(x) cor(x[, y], x[,'Yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c(y, 'Yhat')], x[,'env'], 
                   function(z) cor(z[,y], z[,'Yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(data[, 'gen'], data[,'env'], length)-2)))/
        (sum(1/(tapply(data[, 'gen'], data[, 'env'], length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'Yhat'] - x[, y])^2))
    slope = lapply(reps, function(x) unname(lm(x[,y] ~ x[,'Yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps_M1_CV1,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }else if (cv == 0){
    reps = list()
    for (j in 1:nrept) {
      init[j] = Sys.time()
      set.seed(seed+j)
      means = data[,c('gen','env',y)]
      means$set = NA
      for (i in unique(data$gen)) {
        means[means$gen == i,'set'] = sample(
          1:nfolds, 
          size = dim(means[means$gen == i,])[1],
          replace = dim(means[means$gen == i,])[1] > nfolds
        )
      }
      means$Yhat = NA
      for (k in unique(means$env)) {
        for (i in 1:nfolds) {
          yNA = means[, y]
          yNA[means$env == k] = NA
          
          crossval = BGLR(y = yNA, ETA = ETA, nIter = niter, burnIn = burnin,
                          verbose = F, ...)
          unlink("*.dat")
          
          means$Yhat[means$env == k & means$set == i] = crossval$yHat[means$env == k & means$set == i]
        }
      }
      reps[[j]] = means
      rm(means, crossval, yNA)
      end[j] = Sys.time()
      setTxtProgressBar(pb, j)
      time = round(lubridate::seconds_to_period(sum(end-init)),0)
      est = nrept * (mean(end[end != 0] - init[init != 0])) - time
      remaining = round(lubridate::seconds_to_period(est), 0)
      cat(paste(" // Execution time:", time,
                " // Time you have to take a nap", remaining), "")
    }
    
    corr = lapply(reps, function(x) cor(x[, y], x[,'Yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c(y, 'Yhat')], x[,'env'], 
                   function(z) cor(z[,y], z[,'Yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(data[, 'gen'], data[,'env'], length)-2)))/
        (sum(1/(tapply(data[, 'gen'], data[, 'env'], length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'Yhat'] - x[, y])^2))
    slope = lapply(reps, function(x) unname(lm(x[,y] ~ x[,'Yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps_M1_CV1,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
    
  }else if (cv == 00){
    set.seed(seed)
    sets = split(
      rep(
        1:nfolds, length(unique(data$gen)) * nrept
      )[order(runif(length(unique(data$gen)) * nrept))],
      f = 1:nrept
    )
    reps = list()
    for (j in 1:nrept) {
      init[j] = Sys.time()
      means = data[,c('gen','env',y)]
      means = merge(means, data.frame(
        gen = unique(data$gen),
        set = sets[[j]]
      ), by = gen)
      means$Yhat = NA
      for (k in  unique(means$env)) {
        for (i in 1:nfolds) {
          yNA = means[, y]
          yNA[means$env == k & means$set == i] = NA
          
          cv = BGLR(y = yNA, ETA = ETA, nIter = niter, burnIn = burnin,
                    verbose = F, ...)
          unlink("*.dat")
          
          means$Yhat[means$env == k & means$set == i] = cv$yHat[means$env == k & means$set == i]
        }
      }
      reps[[j]] = means
      rm(means, crossval, yNA)
      end[j] = Sys.time()
      setTxtProgressBar(pb, j)
      time = round(lubridate::seconds_to_period(sum(end-init)),0)
      est = nrept * (mean(end[end != 0] - init[init != 0])) - time
      remaining = round(lubridate::seconds_to_period(est), 0)
      cat(paste(" // Execution time:", time,
                " // Time you have to take a nap", remaining), "")
    }
    
    corr = lapply(reps, function(x) cor(x[, y], x[,'Yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c(y, 'Yhat')], x[,'env'], 
                   function(z) cor(z[,y], z[,'Yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(data[, 'gen'], data[,'env'], length)-2)))/
        (sum(1/(tapply(data[, 'gen'], data[, 'env'], length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'Yhat'] - x[, y])^2))
    slope = lapply(reps, function(x) unname(lm(x[,y] ~ x[,'Yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps_M1_CV1,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }
}






