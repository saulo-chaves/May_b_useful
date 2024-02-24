#' Modified icREML
#' 
#' All credits goes to Professor Verbyla, who created this function

icREML <- function(fm, scale=1) {
  if(!is.list(fm)) stop("Models need to be in a list\n")
  if(is.null(names(fm))) namesfm <- paste("fm", 1:length(fm))
  else namesfm <- names(fm)
  require(asreml)
  asreml.options(Cfixed = TRUE, gammaPar=FALSE)
  fm <- lapply(fm, function(el) {
    if(is.null(el$Cfixed)) {
      out <- update(el, maxit=10) }
    else out <- el
    out})
  logl <- lapply(fm, function(el) el$loglik)
  summ <- lapply(fm, function(el) summary(el, coef=TRUE)$coef.fixed)
  which.X0 <- lapply(summ, function(el) !is.na(el[, 'z.ratio']))
  p.0 <- lapply(which.X0, function(el) sum(el))
  Cfixed <- lapply(fm, function(el) el$Cfixed)
  logdet <- lapply(1:length(fm), function(el, Cfixed, which.X0, scale) {
    log(prod(svd(as.matrix(scale*Cfixed[[el]][
      names(which.X0[[el]][which(which.X0[[el]] != 0)]), 
      names(which.X0[[el]][which(which.X0[[el]] != 0)])
    ]))$d))
  }, Cfixed, which.X0, scale)
  vparam <- lapply(fm, function(el) summary(el)$varcomp)
  q.0 <- lapply(vparam, function(el) sum(!(el$bound == "F" | 
                                             el$bound == "B")))
  b.0 <- lapply(vparam, function(el) sum(el$bound == 'F' | 
                                           el$bound == 'B'))
  logl <- lapply(1:length(fm), function(el, logl, logdet, p.0) {
    logl[[el]] - logdet[[el]]/2}, logl, logdet,p.0)
  aic <- unlist(lapply(1:length(fm), function(el, logl, p.0, q.0) {
    -2*logl[[el]] + 2*(p.0[[el]] + q.0[[el]])}, logl, p.0, q.0))
  bic <- unlist(lapply(1:length(fm), function(el, logl, p.0, q.0, fm) {
    -2*logl[[el]] + log(fm[[el]]$nedf+p.0[[el]])*(p.0[[el]] + q.0[[el]])},
    logl, p.0, q.0, fm))
  results <- data.frame(model=namesfm, loglik = unlist(logl), p=unlist(p.0),
                        q=unlist(q.0), b = unlist(b.0), AIC = aic, 
                        BIC = bic, logdet=unlist(logdet))
  row.names(results) <- 1:dim(results)[1]
  invisible(results)
}