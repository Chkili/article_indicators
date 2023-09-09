S1<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/ENA/S1_RES.ENA.rds")
S2<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/ENA/S2_RES.ENA.rds")
S2<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/ENA/S2_RES.ENA.rds")
S1<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/ENA/S1_RES.ENA.rds")

DH1<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/DH/resultS200000/DH1.rds")
DH2<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/DH/resultS200000/DH2.rds")
DH3<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/DH/resultS200000/DH3.rds")
DH4<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/DH/resultS200000/DH4.rds")

df_S1 <- S1$Sub[,c(1,6,7,15,16)]
df_S1 <- as.data.frame(df_S1)

colnames(S2$Sub)

df_S2 <- S2$Sub[,c(1,6,7,15,16)]
df_S2 <- as.data.frame(df_S2)

r2=NULL

ENAnames <- colnames(df_S2)

for (c in 1 : 5){ # choisir les colonnes que tu veux (les ENA)
  
  x <- df_S1[1:100,c] # modèle de ref du cliff, tu compare donc S2 à S2
  y <- df_S2[1:100,c]
  x.name <- paste0("S2_", ENAnames[c])
  y.name <- paste0("S2_", ENAnames[c])
  results <- as.data.frame(matrix(NA,nrow=1,ncol=8))
  colnames(results) <- c("N","S2","S2","ENA","Cliff estimate","conf.int_min","conf.int_max","magnitude")
  namei<-ENAnames[c]
  cliff.result <- cliff.delta(d=x,f=y)
  results[,1:4] <- c(c,paste(x.name),paste(y.name),namei)
  results[,5:8] <-c(cliff.result$estimate, cliff.result$conf.int,paste(cliff.result$magnitude))
  r2 = rbind(r2, data.frame(results))
  gc()
  assign(paste0("S1_S2_ENA_Cliff"), r2) # la tu créer une dataframe avec tes resultats
}

saveRDS(S1_S2_ENA_Cliff, file="C:/Users/DELL/Documents/Ancien/modelisation gabes/Cliff/S1_S2_ENA_Cliff.rds")
