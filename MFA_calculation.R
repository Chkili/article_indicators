
# calculation of 3000 randomly selected values LIM solutions 

filesr   = list.files("C:/Users/DELL/Documents/Ancien/modelisation gabes/MFA/new_mfa_aleat_env/new_new", pattern = "_variables.rds", full.names = TRUE)


S1<-c(1:300000)
indice<-sample(S1,size=3000,replace=FALSE)
S1<-setdiff(S1,indice)
indice2<-sample(S1,3000)
S1<-setdiff(S1,indice2)
indice3<-sample(S1,3000)
i=1
for (i in 1:length(filesr)){
  print(i)
  fa<-filesr[i] 
  mat=readRDS(fa)
  name = sub("_variables","",tools::file_path_sans_ext(basename(fa)))
  data<-as.data.frame(mat[indice,]) 
 
  
  data<-select(data,  c("gpp->mic","gpp->nan","gpp->pic","doc->bac","mic->los","nan->los","mes->los","R4","R6","R7","R8","TST","AMI","APL","DH","FCI","ACratio","Ninorg","Norg","Pinorg","Porg","SiOH"))%>%
    mutate(eco=str_extract(name,"^[^_]*")) 
  
  if (i>1){
    dataa<-rbind(dataa,data)
  } else{
    dataa<-data
  }
  rm(data)
}


##plot_MFA##

#plot-variables
  
  varr<- dataa
  colnames(varr)
  res.mfa <- MFA(varr, 
                 group = c(3,1,3,4,6,5,1),
                 type = c("s","s","s","s","s","s","n"),
                 name.group = c("primary_prod", "bacprod","Export","Ratio","ENA","Env","web"),
                 num.group.sup = NULL,
                 graph = FALSE)
  
  print(res.mfa)
  
  
  va<- res.mfa$quanti.var   ; va
  vaa<-va$contrib         ; vaa
  #vaaa<-vaa$col.abs     ; vaaa
  #fviz_mfa_var(res.mfa, repel= TRUE)
  X11();plot(res.mfa,choix="var",axes = c(1,2),cex=1,repel=TRUE)
  
  
  #plot_station#
  
  Systeme<-as.factor(c(rep("S1",3000),rep("S2",3000),rep("S3",3000),rep("S4",3000)))
  head(Systeme)
  par(mfrow=c(1,3))
  #s.class(res.mfa$ind$coord[,1:2],fac=Systeme)
  X11();s.class(res.mfa$ind$coord[,c(1,2)], fac=Systeme ,xax = 1, yax = 2,cstar = 1,
                cellipse = 3, axesell = TRUE,  
                cpoint = 1,col=c("#00AFBB", "#E7B800", "#CC79A7","#52854C"), 
                pch = 20)
  group <- get_mfa_var(res.mfa, "quanti.var")
  
  
 


