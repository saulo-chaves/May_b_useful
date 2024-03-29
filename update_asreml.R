#' Author: Saulo Chaves (saulo.chaves at ufv.br)

up.mod = function(model){
  require(asreml)
  repeat{
    if(any(na.exclude(model$vparameters.pc) >= 1)) {
      model = suppressWarnings(update(model))
    }else{
      break
    }
  }
  return(model)
}
