fa.outs = function(model, name.env, name.gen){
  
  require(asreml)
  
  print("1. Starting fa.outs()")
  data = as.data.frame(model$mf)
  num.env = nlevels(data[, name.env])
  num.gen = nlevels(data[, name.gen])
  
  print("2. Extracting summary")
  sum = summary(model)
  
  print("3. Extracting loadings")
  load = sum$varcomp[grep('fa\\d+', rownames(sum$varcomp)),]
  
  print("4. Extracting variances")
  var = sum$varcomp[grep('var', rownames(sum$varcomp)), 1]
  
  print("5. Preparing loading matrix")
  load$fa = regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings = do.call(cbind, lapply(split(load[,c('fa','component')], f = load$fa), 
                                       function(x) x[,'component']))
  rownames(mat.loadings) = levels(data[, name.env])
  
  print("6. Performing SVD on loading matrix")
  svd.lambda = svd(mat.loadings)
  
  print("7. Adjusting SVD signs if necessary")
  if(sum(svd.lambda$u[,1] < 0)/nrow(svd.lambda$u) > 0.5){
    svd.lambda$u = -1 * svd.lambda$u
    svd.lambda$v = -1 * svd.lambda$v
  }
  
  print("8. Creating rotated loadings")
  mat.loadings.star = svd.lambda$u
  colnames(mat.loadings.star) = colnames(mat.loadings)
  rownames(mat.loadings.star) = rownames(mat.loadings)
  
  print("9. Extracting random coefficients (scores)")
  scores = data.frame(summary(model, coef = T)$coef.random)[
    grep('Comp', rownames(data.frame(summary(model, coef = T)$coef.random))),
  ]
  
  print("10. Preparing scores vector")
  scores$fa = regmatches(rownames(scores), regexpr("Comp\\d+", rownames(scores)))
  scor.vec = do.call(c, lapply(split(scores, f = scores$fa), function(x) x[,'solution']))
  
  print("11. Defining D matrix")
  D = diag(svd.lambda$d^2, nrow = length(svd.lambda$d))
  
  print("12. Dimensions before kronecker")
  print(dim(sqrt(D) %*% t(svd.lambda$v)))
  print(num.gen)
  print(length(scor.vec))
  
  print("13. Rotating scores via kronecker product")
  scor.vec.star = kronecker(sqrt(D) %*% t(svd.lambda$v), diag(num.gen)) %*% scor.vec
  
  print("14. Building rotated score matrix")
  scor.mat.star = matrix(
    scor.vec.star,
    nrow = num.gen,
    ncol = length(unique(scores$fa)),
    dimnames = list(levels(data[, name.gen]), unique(load$fa))
  )
  
  print("14.5 Building genetic covariance matrix Gvcov")
  Gvcov = mat.loadings.star %*% tcrossprod(D, mat.loadings.star) + diag(var)
  
  print("15. Building genetic correlation matrix Gcor")
  Gcor = cov2cor(Gvcov)
  
  
  print("16. Calculating total explained variance")
  expvar = sum(diag(mat.loadings.star %*% tcrossprod(D, mat.loadings.star)))/
    sum(diag(Gvcov))
  
  print("17. Calculating per-environment explained variance")
  expvar.j = matrix(nrow = ncol(mat.loadings), ncol = num.env, 
                    dimnames = list(colnames(mat.loadings), levels(data[, name.env])))
  for (i in 1:nrow(expvar.j)) {
    for (j in 1:ncol(expvar.j)) {
      expvar.j[i,j] = 100 * mat.loadings.star[j,i]^2 * diag(D)[i]/
        (sum(mat.loadings.star[j,]^2 * diag(D)) + var[j])
    }
  }
  
  print("18. Calculating average semivariances")
  lambdacross = mat.loadings.star %*% tcrossprod(D, mat.loadings.star)
  fullsv = list()
  for (i in levels(data[,name.env])) {
    semivar = matrix(nrow = num.env, ncol = 3,
                     dimnames = list(levels(data[,name.env]),
                                     c('i', 'j', 'semivar')))
    for (j in levels(data[,name.env])) {
      semivar[j,] = c(i, j, .5 * (lambdacross[i,i] + lambdacross[j,j]) -
                        lambdacross[i,j])
    }
    fullsv[[i]] = semivar
    rm(semivar)
  }
  
  print("19. Summarizing semivariances for ASV")
  fullsv = data.frame(do.call(rbind, fullsv))
  fullsv$semivar = as.numeric(fullsv$semivar)
  fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[,-1]
  colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv)- 1)) * sum(fullsv[upper.tri(fullsv)])
  rm(fullsv)
  
  print("20. Calculating G-based semivariances for ASV")
  fullsv = list()
  for (i in levels(data[,name.env])) {
    semivar = matrix(nrow = num.env, ncol = 3,
                     dimnames = list(levels(data[,name.env]),
                                     c('i', 'j', 'semivar')))
    for (j in levels(data[,name.env])) {
      semivar[j,] = c(i,j,
                      0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
    }
    fullsv[[i]] = semivar
  }
  
  print("21. Summarizing G-based semivariances")
  fullsv = data.frame(do.call(rbind, fullsv))
  fullsv$semivar = as.numeric(fullsv$semivar)
  fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[,-1]
  colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  G_ASV = 2/(nrow(fullsv) * (nrow(fullsv)- 1)) * sum(fullsv[upper.tri(fullsv)])
  ASVR = lambda_ASV/G_ASV
  
  print("22. Predicting marginal and conditional BLUPs")
  modpred = predict(model, classify = paste(name.gen, name.env, sep = ':'), sed = T,pworkspace=300e6)
  blups = modpred$pvals
  blups = blups[,-5]
  temp = data.frame(
    name.env = rep(levels(data[, name.env]), each = nlevels(data[, name.gen])),
    name.gen = rep(levels(data[, name.gen]), times = nlevels(data[, name.env])), 
    marginal = kronecker(mat.loadings.star, diag(num.gen)) %*% scor.vec.star
  )
  colnames(temp) = c(name.env, name.gen, 'marginal')
  blups = merge(blups, temp, by = c(name.env, name.gen))
  colnames(blups)[which(colnames(blups) == 'predicted.value')] = 'conditional'
  
  print("23. Calculating environment-wise generalized heritabilities")
  colnames(modpred$sed) = rownames(modpred$sed) = paste(modpred$pvals[,1], modpred$pvals[,2], sep = '_')
  
  H2 = NULL
  for (i in levels(data[, name.env])) {
    vd = (modpred$sed[grep(i, rownames(modpred$sed)), grep(i, colnames(modpred$sed))])^2
    H2[i] = 1-(mean(vd[upper.tri(vd ,diag = F)])/(2*diag(Gvcov)[i]))
  }
  
  print("24. Wrapping results into list")
  results = list('rot.loads' = mat.loadings.star, 
                 'Gvcov' = Gvcov, 
                 'Gcor' = Gcor, 
                 'diagnostics' = c(
                   expvar = round(expvar*100, 3),
                   ASVR = round(ASVR * 100, 3),
                   aic = round(sum$aic, 3), 
                   bic = round(sum$bic, 3),
                   logl = round(sum$loglik, 3)
                 ),
                 'expvar_j' = expvar.j,
                 "rot.scores" = scor.mat.star, 
                 'blups' = blups,
                 'H2' = H2)
  
  print("25. Returning results - done!")
  return(results)
}
