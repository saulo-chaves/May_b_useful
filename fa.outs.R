fa.outs = function(model, name.env, name.gen){
  
  require(asreml)
  
  data = as.data.frame(model$mf)
  num.env = nlevels(data[, name.env])
  num.gen = nlevels(data[, name.gen])
  sum = summary(model)
  load = sum$varcomp[grep('fa\\d+', rownames(sum$varcomp)),]
  var = sum$varcomp[grep('var', rownames(sum$varcomp)),1]
  
  
  load$fa = regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings = do.call(cbind, lapply(split(load[,c('fa','component')], f = load$fa), 
                                       function(x) x[,'component']))
  rownames(mat.loadings) = levels(data[, name.env])
  mat.loadings.star = mat.loadings %*% svd(mat.loadings)$v
  if(sum(mat.loadings.star[,1] < 0)/dim(mat.loadings.star)[1] > .2) mat.loadings.star = mat.loadings.star * -1
  colnames(mat.loadings.star) = unique(load$fa)
  
  Gvcov = tcrossprod(mat.loadings.star) + diag(var)
  Gcor = cov2cor(Gvcov)
  
  
  expvar = sum(diag(tcrossprod(mat.loadings.star)))/sum(diag(Gvcov))
  expvar.j = matrix(nrow = ncol(mat.loadings), ncol = num.env, 
                    dimnames = list(colnames(mat.loadings), levels(data[, name.env])))
  for (i in 1:nrow(expvar.j)) {
    for (j in 1:ncol(expvar.j)) {
      expvar.j[i,j] = 100 * mat.loadings.star[j,i]^2/(sum(mat.loadings.star[j,]^2) + var[j])
    }
  }
  
  
  scores = data.frame(summary(model, coef = T)$coef.random)[
    grep('Comp', rownames(data.frame(summary(model, coef = T)$coef.random))),
  ]
  scores$fa = regmatches(rownames(scores), regexpr("Comp\\d+", rownames(scores)))
  scor.vec = do.call(c, lapply(split(scores, f = scores$fa), function(x) x[,'solution']))
  scor.vec.star = -kronecker(t(svd(mat.loadings)$v), diag(num.gen)) %*% scor.vec
  
  scor.mat.star = matrix(scor.vec.star, nrow = num.gen, 
                         ncol = length(unique(scores$fa)),
                         dimnames = list(levels(data[, name.gen]), unique(load$fa)))
  
  modpred = predict(model, classify = paste(name.env, name.gen, sep = ':'))
  blups = modpred$pvals
  blups = blups[,-5]
  blups$marginal = kronecker(mat.loadings.star, diag(num.gen)) %*% scor.vec.star
  colnames(blups)[which(colnames(blups) == 'predicted.value')] = 'conditional'
  
  
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
                 'blups' = blups)
  
  return(results)
}










