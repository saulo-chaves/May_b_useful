##' Function crossvalid
##'
##' @title Cross-validation schemes
##' 
##' @description
##' This functions performs parallelized cross-validations for obtaining the predictive ability
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
##' @import BGLR parallel
##' 
##' @author Saulo F. S. Chaves (saulo.chaves at ufv.br)
##'

crossvalid = function(cvdata,  ETA, cv = 1, results = NULL,...){
  
  cvdata = cvdata

  if (cv == 1){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo['nrept']), 
                             f = 1:ceiling(cvdata$cvinfo['nrept']/ parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo['nrept']), 
              f = 1:ceiling(cvdata$cvinfo['nrept']/parallel::detectCores(logical = F)))[[i]]
        ], function(z){
        lapply(1:cvdata$cvinfo['nfolds'], function(x){
          
          yNA = z$trait
          yNA[z$set == x] = NA
          
          crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo['niter'], 
                                burnIn = cvdata$cvinfo['burnin'], verbose = F)
          unlink("*.dat")
          
          cbind(
            z[z$set == x,], 
            yhat = crossval$yHat[z$set == x]
          )
        })
      })
    }
    
    parallel::stopCluster(cl)
    
    reps = unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    )
    
    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(cvdata$cvdata[[1]]$gen, 
                               cvdata$cvdata[[1]]$env, length)-2)))/
        (sum(1/(tapply(cvdata$cvdata[[1]]$gen, 
                       cvdata$cvdata[[1]]$env, length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }else if (cv == 2){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo['nrept']), 
                             f = 1:ceiling(cvdata$cvinfo['nrept']/parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo['nrept']), 
              f = 1:ceiling(cvdata$cvinfo['nrept']/parallel::detectCores(logical = F)))[[i]]
        ], function(z){
        lapply(1:cvdata$cvinfo['nfolds'], function(x){
          
          yNA = z$trait
          yNA[z$set == x] = NA
          
          crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo['niter'],
                                burnIn = cvdata$cvinfo['burnin'], verbose = F, ...)
          unlink("*.dat")
          
          cbind(
            z[z$set == x,], 
            yhat = crossval$yHat[z$set == x]
          )
        })
      })
    }

    parallel::stopCluster(cl)
    
    reps = unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    )

    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(cvdata$cvdata[[1]]$gen, 
                               cvdata$cvdata[[1]]$env, length)-2)))/
        (sum(1/(tapply(cvdata$cvdata[[1]]$gen, 
                       cvdata$cvdata[[1]]$env, length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }else if (cv == 0){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo['nrept']),
                             f = 1:ceiling(cvdata$cvinfo['nrept']/parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, X = cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo['nrept']),
              f = 1:ceiling(cvdata$cvinfo['nrept']/parallel::detectCores(logical = F)))[[i]]
        ], fun = function(z){
        lapply(unique(cvdata$cvdata[[1]]$env), function(w){
          lapply(1:cvdata$cvinfo['nfolds'], function(x){
            
            yNA = z$trait
            yNA[z$env == w] = NA
            
            crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo['niter'],
                                  burnIn = cvdata$cvinfo['burnin'], verbose = F, ...)
            unlink("*.dat")
            
            cbind(
              z[z$set == x & z$env == w,], 
              yhat = crossval$yHat[z$set == x & z$env == w]
            )
          })
        })
      })
    }
    
    parallel::stopCluster(cl)
    
    reps = lapply(unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    ), function(y){do.call(rbind,y)})
    
    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(cvdata$cvdata[[1]]$gen, 
                               cvdata$cvdata[[1]]$env, length)-2)))/
        (sum(1/(tapply(cvdata$cvdata[[1]]$gen, 
                       cvdata$cvdata[[1]]$env, length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
    
  }else if (cv == 00){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo['nrept']), 
                             f = 1:ceiling(cvdata$cvinfo['nrept']/parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, X = cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo['nrept']),
              f = 1:ceiling(cvdata$cvinfo['nrept']/parallel::detectCores(logical = F)))[[i]]
        ], fun = function(z){
        lapply(unique(cvdata$cvdata[[1]]$env), function(w){
          lapply(1:cvdata$cvinfo['nfolds'], function(x){
            
            yNA = z$trait
            yNA[z$env == w & z$set == x] = NA
            
            crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo['niter'],
                                  burnIn = cvdata$cvinfo['burnin'], verbose = F, ...)
            unlink("*.dat")
            
            cbind(
              z[z$set == x & z$env == w,], 
              yhat = crossval$yHat[z$set == x & z$env == w]
            )
          })
        })
      })
    }
    
    parallel::stopCluster(cl)
    
    reps = lapply(unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    ), function(y){do.call(rbind,y)})
 
    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    wapa = lapply(corr_env, function(x){
      (sum((x/(1-x^2))/(tapply(cvdata$cvdata[[1]]$gen, 
                               cvdata$cvdata[[1]]$env, length)-2)))/
        (sum(1/(tapply(cvdata$cvdata[[1]]$gen,
                       cvdata$cvdata[[1]]$env, length)-2)))
    })
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
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



### Add description

assignation = function(data, y, gen, env, seed, nfolds, nrept, cv = 1, 
                       niter, burnin){
  
  data = data
  data$gen = data[,gen]
  data$env = data[,env]
  data$trait = data[,y]
  
  if (cv == 1){
    set.seed(seed)
    sets = split(
      rep(
        1:nfolds, length(unique(data$gen)) * nrept
      )[order(runif(length(unique(data$gen)) * nrept))],
      f = 1:nrept
    )
    cvdata = lapply(sets, function(x){
      cvdata = data[,c('gen','env','trait')]
      cvdata = merge(cvdata, data.frame(
        gen = unique(data$gen),
        set = x
      ), by = gen)
    })
    names(cvdata) = paste0('R',1:nrept)
    
  }else if (cv == 2){
    cvdata = list()
    for (j in 1:nrept) {
      set.seed(seed+j)
      cvdata[[j]] = data[,c('gen','env','trait')]
      cvdata[[j]]$set = NA
      for (i in unique(data$gen)) {
        cvdata[[j]][cvdata[[j]]$gen == i,'set'] = sample(
          1:nfolds, 
          size = dim(cvdata[[j]][cvdata[[j]]$gen == i,])[1],
          replace = dim(cvdata[[j]][cvdata[[j]]$gen == i,])[1] > nfolds
        )
      }
    }
    names(cvdata) = paste0('R',1:nrept)
  }
  cvinfo = c(nfolds = nfolds, nrept = nrept, cv = cv, niter = niter, burnin = burnin)
  
  return(list(cvdata = cvdata, cvinfo = cvinfo))
}


