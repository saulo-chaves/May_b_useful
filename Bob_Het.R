library(data.table)
library(AGHmatrix)
library(tidyverse)


data_SNP<- fread("02_output\\tables\\SNP_DATA.txt")

data_SNP <- data_SNP %>% remove_rownames %>% column_to_rownames(var="V1")

data_SNP[1:5, 1:5]


data_SNP[data_SNP== 0] <- "--"
data_SNP[data_SNP== '00'] <- "--"
a <- data_SNP[1:500, 1:500]


dim(data_SNP)

bob_het = function(data, NA_lim, HET_lim){   # O SNP que passar do limite de dados perdidos "NA_lim" ou do limete de heterozigotos "HET_lim, será retirado do banco
  if(!require(dplyr)){install.packages("dplyr")}
  
  print(paste("Started at", strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S")))
  
  sample = data_SNP
  
  # Quantidade de cada combinação para um dado SNP
  sample$nHET = rowSums(sample != "CC" & sample != "AA" &        # criando colunas no df
                          sample != "TT" & sample!= "GG" & sample!= "--") # Tudo que for diferente de homozigoto fica nesta coluna
  sample$nAA = rowSums(sample == "AA")  # Tudo que for homozigoto fica nestas coluna
  sample$nTT = rowSums(sample == "TT")
  sample$nGG = rowSums(sample == "GG")
  sample$nCC = rowSums(sample == "CC")
  sample$nNA = rowSums(sample == "--")
  
  abc <- sample %>% select( nAA, nTT, nGG, nCC)
  

  abc$MAX = apply(abc, 1, which.max) 
  max2<- abc$MAX
  #max2[max2 == 1] <- "HET"
  max2[max2 == 1] <- "AA"
  max2[max2 == 2] <- "TT"
  max2[max2 == 3] <- "GG"
  max2[max2 == 4] <- "CC"
  

  # A combinação mais frequente de cada linha
  
  # The sapply in R is a built-in function that applies a function to all the input elements.
  
  # Diagnósticos
  print(paste("Data has", length(which(sample$nAA == ncol(data) | 
                                         sample$nTT == ncol(data) | 
                                         sample$nGG == ncol(data) |
                                         sample$nCC == ncol(data) | 
                                         sample$nNA == ncol(data) | 
                                         sample$nHET == ncol(data))),"monomorphic markers"))
  
  print(paste("Data has", length(which(sample$nNA/ncol(data)*100 > NA_lim)), 
              "markers with more than", NA_lim,"% of missing values"))
  
  print(paste("Data has", length(which(sample$nHET/ncol(data)*100 > HET_lim)), 
              "markers with more than", HET_lim, "% of heterozigotic loci"))
  
  
  # Filtros (Sample filtered)- O filtro deixa quem cumpre a condição
  # Ex: se o nAA for diferente do número total de colunas, significa que aquele SNP não é monomórfico para AA
  
  SF = sample %>%  
    filter(nAA != ncol(data) & nTT != ncol(data) & nGG != ncol(data) &
             nCC != ncol(data) & nNA != ncol(data) & nHET != ncol(data)) %>% 
    filter((nNA/ncol(data))*100 <= NA_lim) %>% filter((nHET/ncol(data))*100 <= HET_lim) %>%  #NA_lim e HET_lim
    select(-c(nAA, nTT, nGG, nCC, nNA, nHET))            # Selecione todos que cumprem as condições (tira os monomórfico, Call Rate e Heterozigotos(quando necessário)
  
 # prev = SF$prev   # Combinações alélicas que prevalecem em cada marca
  # SF = select(SF, -prev)  # Matriz de marcadores não codificada, más já filtrada
  
  
  print(paste('Based on the provided limits,',
              nrow(sample) - nrow(SF), "markers were discarded, so we are left with",
              nrow(SF), "markers"))
  SF = as.matrix(SF) # tranforma o SF em matrix
  
  SNP = SF
  
  # Codificando a matriz
  
  for(i in 1:nrow(SF)){
    
    SNP[i,][SNP[i,] != "--" & SNP[i,] != "AA" & 
              SNP[i,] != "TT" & SNP[i,] != "CC" & 
              SNP[i,] != "GG"]  = 1
    SNP[i,][SNP[i,] == "--"] = -9
    SNP[i, ] = ifelse(SNP[i, ] == max2[i], 2, SNP[i, ])  ####
    
    SNP[i,][SNP[i,] != "-9" & SNP[i,] != "1" & SNP[i,] != "2"] = 0
    
  }
  
  print(paste("Ended at", strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S")))
  return(SNP)
  
}

teste = bob_het(data_SNP, NA_lim = 5, HET_lim = 100)  
dim(teste)

# No NA_lim coloca-se o limite de dados perdidos aceitável para um dado SNP, 
# por exemplo se O SNP tiver mais do que 5% de dados perdido, excluá-o.

# No HET_lim coloca-se o limite de heterozigotos que é aceitável para um dado SNP, 
# por exemplo se o limite for 100, nenhum SNP será excluído. Usa-se HET_lim igual a 0, 
# quando os SNPs são oriundo de linhagens puras.

# Retirar genótipos

# b <- as.matrix(colSums(teste == -9))
# colnames(b) = "b"
# dim(b)
# 
# p_corte = nrow(teste) * 0.05
# p_corte
# 
# c = as.data.frame(b) %>% filter(b < p_corte)
# c <- t(c)
# 
# teste.indfilt = teste[,intersect(colnames(c),colnames(teste))]
# 
# dim(teste.indfilt)


# MAF ---------------------------------------------------------------------

MAF = NULL
p = NULL
q = NULL
for(j in 1:nrow(teste)){
  Z = length(which("-9" == teste[j,]))
  D = length(which(2 == teste[j,]))/(ncol(teste)-Z)
  H = length(which(1 == teste[j,]))/(ncol(teste)-Z)
  R = length(which(0 == teste[j,]))/(ncol(teste)-Z)
  p[j] = D + (H/2)
  q[j] = R + (H/2)
  MAF[j] = min(p[j],q[j])
  
  print(j)
}

teste_df = as.data.frame(teste)
teste_df$MAF = MAF

teste_df = teste_df %>% filter(MAF >= 0.05)   # Pode-se testar vários valores de MAF
teste_df = teste_df %>%  select(-MAF)
dim(teste_df)

write.table(teste_df, "02_output//tables//FILT_GENO_DATA.txt", sep= " ")

SNPmat = read.table("02_output//tables//FILT_GENO_DATA.txt")

SNPmat[1:10, 1:10]

SNPmat = t(as.matrix(SNPmat))

A_mat = Gmatrix(SNPmat, method = "VanRaden")
heatmap(A_mat)

write.table(A_mat, "02_output//tables//A_MAT_GRANDIS.txt", sep = " ")

D_mat = Gmatrix(SNPmat, method = "Vitezica")
heatmap(D_mat)

write.table(D_mat, "02_output//tables//D_MAT_GRANDIS.txt", sep = " ")



