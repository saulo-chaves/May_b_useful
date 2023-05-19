##' Function fa.outs
##'
##' @title Extracting outputs from factor analytic mixed models
##' 
##' @description
##' This function compiles the outputs of factor analytic mixed models fitted using 
##' [asreml::asreml()], and use them to compute several useful parameters, such as the 
##' genetic (co)variances, genetic correlations between environments, environment-wise 
##' generalized heritabilities and marginal eBLUPs.
##' 
##' @param model An object of class `asreml` containing the results of the fitted 
##' factor analytic mixed model.
##' @param name.env A string that represents the name of the "environment" factor of the model
##' @param name.gen A string that represents the name of the "genotype" factor of the model
##' 
##' @return The function returns a list with: 
##' \itemize{
##' \item \code{rot.loads}: a matrix with the rotated loadings 
##' \item \code{Gvcov}: the full genetic variance-covariance matrix 
##' \item \code{Gcor}: the matrix of genetic correlations between environments 
##' \item \code{diagnostics}: goodness-of-fit diagnostics: percentage of explained variance, 
##' AIC, BIC and log(L).
##' \item \code{expvar.j}: a matrix containing the percentage of explained variance per environment, 
##' by each factor.
##' \item \code{rot.scores}: a matrix with the rotated scores.
##' \item \code{blups}: a dataframe with the conditional eBLUPs (and their standard error), and 
##' marginal eBLUPs (i.e., disregarding the specific effects).
##' \item \code{H2}: a vector with the environment-wise generalized heritabilities.
##' }
##' 
##' @details
##' 
##' If this turn into a package someday, I will add some details here.
##' 
##' 
##' @author Saulo F. S. Chaves (saulo.chaves at ufv.br)
##'

fa.outs = function(model, name.env, name.gen){
  
  require(asreml)
  
  data = as.data.frame(model$mf)
  num.env = nlevels(data[, name.env])
  num.gen = nlevels(data[, name.gen])
  sum = summary(model)
  load = sum$varcomp[grep('fa\\d+', rownames(sum$varcomp)),]
  var = sum$varcomp[grep('var', rownames(sum$varcomp)),1]
  
  # Loadings rotation via SVD
  load$fa = regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings = do.call(cbind, lapply(split(load[,c('fa','component')], f = load$fa), 
                                       function(x) x[,'component']))
  rownames(mat.loadings) = levels(data[, name.env])
  mat.loadings.star = mat.loadings %*% svd(mat.loadings)$v
  if(sum(mat.loadings.star[,1] < 0)/dim(mat.loadings.star)[1] > .2) mat.loadings.star = mat.loadings.star * -1
  
  # Full variance-covariance genetic matrix and genetic correlation matrix
  Gvcov = tcrossprod(mat.loadings.star) + diag(var)
  Gcor = cov2cor(Gvcov)
  
  # Percentage of explained variance
  expvar = sum(diag(tcrossprod(mat.loadings.star)))/sum(diag(Gvcov))
  expvar.j = matrix(nrow = ncol(mat.loadings), ncol = num.env, 
                    dimnames = list(colnames(mat.loadings), levels(data[, name.env])))
  for (i in 1:nrow(expvar.j)) {
    for (j in 1:ncol(expvar.j)) {
      expvar.j[i,j] = 100 * mat.loadings.star[j,i]^2/(sum(mat.loadings.star[j,]^2) + var[j])
    }
  }
  
  # Scores rotation via SVD
  scores = data.frame(summary(model, coef = T)$coef.random)[
    grep('Comp', rownames(data.frame(summary(model, coef = T)$coef.random))),
  ]
  scores$fa = regmatches(rownames(scores), regexpr("Comp\\d+", rownames(scores)))
  scor.vec = do.call(c, lapply(split(scores, f = scores$fa), function(x) x[,'solution']))
  scor.vec.star = kronecker(t(svd(mat.loadings)$v), diag(num.gen)) %*% scor.vec
  if(sum((mat.loadings %*% svd(mat.loadings)$v)[,1] < 0)/dim((mat.loadings %*% svd(mat.loadings)$v))[1] > .2) scor.vec.star = scor.vec.star * -1
  scor.mat.star = matrix(scor.vec.star, nrow = num.gen, ncol = length(unique(scores$fa)),
                         dimnames = list(levels(data[, name.gen]), unique(load$fa)))
  
  # eBLUPs
  modpred = predict(model, classify = paste(name.env, name.gen, sep = ':'), sed = T)
  blups = modpred$pvals
  blups = blups[,-5]
  blups$marginal = kronecker(mat.loadings.star, diag(num.gen)) %*% scor.vec.star
  colnames(blups)[which(colnames(blups) == 'predicted.value')] = 'conditional'
  
  # Environment-wise generalized heritabilities
  colnames(modpred$sed) = rownames(modpred$sed) = paste(modpred$pvals[,1], modpred$pvals[,2], sep = '_')
  
  H2 = NULL
  for (i in levels(data[, name.env])) {
    vd = (modpred$sed[grep(i, rownames(modpred$sed)), grep(i, colnames(modpred$sed))])^2
    H2[i] = 1-(mean(vd[upper.tri(vd ,diag = F)])/(2*diag(Gvcov)[i]))
  }
  
  results = list('rot.loads' = mat.loadings.star, 
                 'Gvcov' = Gvcov, 
                 'Gcor' = Gcor, 
                 'diagnostics' = c(
                   expvar = round(expvar*100,3), 
                   aic = round(sum$aic,3), 
                   bic = round(sum$bic,3),
                   logl = round(sum$loglik,3)
                 ),
                 'expvar_j' = expvar.j,
                 "rot.scores" = scor.mat.star, 
                 'blups' = blups,
                 'H2' = H2)
  
  return(results)
}

