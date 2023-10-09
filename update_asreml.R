up.mod = function(model){
  repeat{
    if(any(na.exclude(model$vparameters.pc) >= 1)) {
      model = update(model)
    }else{
      break
    }
  }
  return(mod)
}