### We built this function for coding and filtering the SNP matrices.
### 0 is for the less frequent allele, 2 is for the most frequent, 1 is for the heterozigote and -9 is for the missing data
### In this version, we assume that your SNPs came from inbred lines
### Your dataframe must have the SNPs in the row and the genotypes in the column
### We are filtering for call rate (missing data > NA_lim is excluded) and frequency of heterozigote (%Heterozigote > HET_lim is excluded)
### MAF filtering will be added to the function

bob_het = function(data, NA_lim, HET_lim){
  if(!require(dplyr)){install.packages("dplyr")}
  
  print(paste("Started at", strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S")))
  
  sample = data
  
  sample$nHET = rowSums(sample != "CC" & sample != "AA" &
                          sample != "TT" & sample!= "GG" & sample!= "--")
  sample$nAA = rowSums(sample == "AA")
  sample$nTT = rowSums(sample == "TT")
  sample$nGG = rowSums(sample == "GG")
  sample$nCC = rowSums(sample == "CC")
  sample$nNA = rowSums(sample == "--")
  sample$prev = names(sapply(apply(sample, 1, table), which.max))
  
  print(paste("Data has", length(which(sample$nAA == ncol(data) | 
                                 sample$nTT == ncol(data) | 
                                 sample$nGG == ncol(data) |
                                 sample$nCC == ncol(data) | 
                                 sample$nNA == ncol(data) | 
                                 sample$nHET == ncol(data))),"monomorphic markers"))
  
  print(paste("Data has", length(which(sample$nNA/ncol(data)*100 >= NA_lim)), 
      "markers with more than", NA_lim,"% of missing values"))
  
  print(paste("Data has", length(which(sample$nHET/ncol(data)*100 > HET_lim)), 
      "markers with more than", HET_lim, "% of heterozigotic loci"))
  
  SF = sample %>%  
    filter(nAA != ncol(data) & nTT != ncol(data) & nGG != ncol(data) &
             nCC != ncol(data) & nNA != ncol(data) & nHET != ncol(data)) %>% 
    filter((nNA/ncol(data))*100 <= NA_lim) %>% filter((nHET/ncol(data))*100 <= HET_lim) %>% 
    select(-c(nAA, nTT, nGG, nCC, nNA, nHET))
  
  prev = SF$prev
  SF = select(SF, -prev)
  print(paste('Based on the provided limits,',
      nrow(sample) - nrow(SF), "markers were discarded, so we are left with",
      nrow(SF), "markers"))
  SF = as.matrix(SF)
  
  SNP = SF
  
  for(i in 1:nrow(SF)){
    
    SNP[i,][SNP[i,] != "--" & SNP[i,] != "AA" & 
              SNP[i,] != "TT" & SNP[i,] != "CC" & 
              SNP[i,] != "GG"]  = 1
    SNP[i,][SNP[i,] == "--"] = -9
    SNP[i, ] = ifelse(SNP[i, ] == prev[i], 2, SNP[i, ])
    
    SNP[i,][SNP[i,] != "-9" & SNP[i,] != "1" & SNP[i,] != "2"] = 0
    
  }

  print(paste("Ended at", strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S")))
  return(SNP)
}


