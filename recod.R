recod = function(data, name.fact, cod){
  stopifnot("cod must be a string" = is.character(cod))
  stopifnot("name.fact must be a string" = is.character(name.fact))
  
  num = length(unique(data[,name.fact]))
  data$ord = 1:nrow(data)
  
  a = data.frame(
    "JJJJ" = unique(data[,name.fact]),
    'XYWFX' = paste0(cod,sprintf(paste0('%0', nchar(num),'d'), seq(1:num)))
  ) 
  
  colnames(a)[which(colnames(a) == 'JJJJ')] = name.fact
  
  data = merge(data, a, by = name.fact, sort=F)
  colnames(data)[which(colnames(data) == "XYWFX")] = paste("cod",name.fact,sep = '.')
  
  data = data[order(data$ord),]
  data = data[,-which(colnames(data)=='ord')]
  rownames(data) = NULL
  
  return(data)
  
}
