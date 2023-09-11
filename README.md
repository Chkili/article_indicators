# Journal of Marine Systems 
This R script comes with the following research article:
It includes (i) a procedure for calculation of the unknown fluxes using the LIM-MCMC method based on the mirror technique defined by Meersche et al. (2009) estimates each unknown flux, (ii) a procedure for estimating ENA indicators and typology ratios and (iv) procedure to create a Multiple Factor Analysis (MFA) Ordination Diagram showing the relationships between ecological indicators, environmental variables and carbon fluxes.

## Data_matrix
[Station1-S1](https://github.com/Chkili/article_indicators/blob/main/S1_mod.txt)
<br>
[Station2-S2](https://github.com/Chkili/article_indicators/blob/main/S2_mod.txt)
<br>
[Station3-S3](https://github.com/Chkili/article_indicators/blob/main/S3_mod.txt)
<br>
[Station4-S4](https://github.com/Chkili/article_indicators/blob/main/S4_mod.txt)

## Load required libraries
```R
library(LIM)
library(limSolve)
library(diagram)
library (shape)
library(MASS)
library(ade4)
library(vegan)
library(NetIndices)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(qgraph)
library(igraph)
library(ggplot2)
library(effsize)
library(stringr)
```

## Calculation of the unknown fluxes

In our study, four food web typology ratios of Sakka Hlaili et al. (2014)  were calculated from the flux data yielded by the models to describe the different interactions between compartments and identify the type of trophic pathway.
```R
 #to set which directory you want to work in: source path
setwd("C:/Users/...") 

#to check the right path
getwd() 

#download declaration file and transformation into ABGH matrix
S1.lim <- Setup("Golfe_S1.txt") 

#pars it's to calculate the deterministic (least squares) (parsimonious) solution (a single value for a flow instead of a set of solutions)
pars <- Ldei(S1.lim) 

#xranges it's to finds the possible ranges ([min,max]) for each unknown.
webranges<- Xranges(S1.lim) 

# calculation of random solution jump=10 iter=300000
S1.X<-Xsample(S1.lim,iter=300000,jmp=10) 
```
## Calculation of Food web typology ratios 
``` R
#typology ratios have been calculated from certain Xsample flows


S1<- S1.X #x.sample

phytTOpro<- S1[,7]+S1[,13]+S1[,18]
D1<- phytTOpro+ S1[,29]

GPPphyt<-S1[,1]+S1[,2]+S1[,3]
resphyt<-S1[,4]+S1[,10]+S1[,16]

Pnet_phyt<- GPPphyt-resphyt
Pnet_bac<-S1[,30]- S1[,27]
Pnet_doc<- S1[,5]+S1[,11]+S1[,17]+S1[,20]+S1[,24]+S1[,28]+S1[,32]
Pnet_det<-S1[,6]+S1[,12]+S1[,21]+S1[,25]
D2<- Pnet_phyt+Pnet_bac+Pnet_doc+Pnet_det
R4<- Pnet_phyt/D2           
summary(R4)


#R6=(Pnetdet + Pnetdoc) /(Pnetpht + Pnetbac + Pnetdet + Pnetdoc)
R6<-(Pnet_det + Pnet_doc)/D2
summary(R6)

#R7=Pnetpico/Pnetpht
Pnet_pico<- S1[,3]-S1[,16]
R7<-Pnet_pico/Pnet_phyt
summary(R7)

#R8=phtTOpro/ (phtTOpro + phtTOmet)
phytTOmet<- S1[,8]+S1[,14]
R8<- phytTOpro/(phytTOpro+phytTOmet)
summary(R8)
########################
#summary ratio##########
########################
summary(R4)
summary(R6)
summary(R7)
summary(R8)
```
### Plot Food web typology ratios 
```R
ratio.S1<-readRDS("C:/Users/.../S1_ratio.rds")
ratio.S2<-readRDS("C:/Users/.../S2_ratio.rds")
ratio.S3<-readRDS("C:/Users/.../S3_ratio.rds")
ratio.S4<-readRDS("C:/Users/.../S4_ratio.rds")

summary(ratio.S1)
summary(ratio.S2)
summary(ratio.S3)
summary(ratio.S4)


S1<-gather(ratio.S1[,c(1:6)], ratios)%>%
  mutate( web="S1")
summary(S1)
S2<-gather(ratio.S2[,c(1:6)], ratios)%>%
  mutate( web="S2")
S3<-gather(ratio.S3[,c(1:6)], ratios)%>%
  mutate( web="S3")
S4<-gather(ratio.S4[,c(1:6)], ratios)%>%
  mutate( web="S4")
#d18<-gather(dd18[,c(1:6)], ratios)%>%
#mutate( web="d18")

ratio<-bind_rows(S1,S2,S3,S4)
ratio
summary(ratio)

X<-filter(ratio,ratios %in% c("R4","R6", "R7", "R8"))%>%
  mutate(ratioss=as.factor(ratios))

summary(X)
p<-ggplot(data=X,aes(x=web,y=value))
x11();p+geom_boxplot()+theme_classic()+facet_wrap(~ratios,scale="free")+
  theme(strip.text.x=element_text(size=16), axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), axis.title =element_text(size=16))+
  labs(x="Web",y="values")
# Prepare a vector of colors with specific color for Nairobi and Eskimo
myColors <- ifelse(levels(X$web)=="d00" , rgb(0.1,0.1,0.7,0.5) , "grey90" )
```

![Typology_ratios](https://github.com/Chkili/article_indicators/blob/main/ratios..jpg)

Spatial variation of the food web typology ratios. 
R4 = Pnetpht/D2 (A), R6 = (PnetDET + PnetDOC)/D2 (B), R7 = PnetPIC/Pnetpht (C) and R8 = pht-PRO/(pht-PRO + pht-MET) (D)
Pnet, net production; pht, total phytoplankton; D2, Pnetpht + PnetBAC + PnetDET + PnetDOC; pht-PRO, consumption of pht by PRO; pht-MET, consumption of pht by MET.

## Calculation of ENA indices

After solution, ENA indices were estimated for each station (i.e. S1, S2, S3 and S4), using functions from
package NetIndices (Soetaert and Kones 2008; Kones, Soetaert, van Oevelen, and Owino 2009). 
Several indices have been calculated, but we have chosen 5 which are : AMI, A/C, TST; APL, FCI
```R
Indices <-  function(S1.lim,    # web specifications
                     S1.X, # random samples
                     ext,        # external compartments (import and export)
                     dead,       # dead compartments (TL=1)
                     specs)      # species for which trophic analysis is wanted
{
  # binary flow matrix
  flowmat    <- S1.lim$Flowmatrix
  flowmatrix <- flowmat
  ii         <- which(flowmat > 0, arr.ind = TRUE) #search all feeds >0
  
  # Random sample
  X          <- S1.X # corresponds to the random solutions
  Indices    <- NULL  # will contain all network indices
  indices    <- NULL  # a subset
  
  for (i in 1:nrow(X))
  {
    
    # generate required flowmatrix for this solution
    flowmatrix[ii] <- X[i,flowmat[ii]]
    Ing<-rowSums   (flowmatrix,dims=1)
    
    # calculate all network indices
    UU<-UncInd    (flowmatrix,Import=ext,Export=ext)
    EF<-EffInd    (flowmatrix,Import=ext,Export=ext)
    AA<-AscInd    (flowmatrix,Import=ext,Export=ext)
    EE<-EnvInd    (flowmatrix,Import=ext,Export=ext)
    GG<-GenInd    (flowmatrix,Import=ext,Export=ext)
    PP<-PathInd   (flowmatrix,Import=ext,Export=ext)
    TT<-TrophInd  (flowmatrix,Import=ext,Export=ext,Dead=dead)
    
    # select the indices that will be investigated
    Ind<-unlist(c(UU,AA[1,],AA[2,],EE, EF, GG, PP,TT))
    
    ind <- c(Ind["T.."],Ind["Cbar"],Ind["Ascendency"],Ind["Overhead"],
             Ind["Capacity"],Ind["ACratio"],Ind["AMI"],Ind["HR"],
             Ind["DR"],Ind["CZ"],Ind["FZ"],Ind["NZ"],Ind["RZ"],
             Ind["ID"],Ind["FCI"],Ind["APL"],Ind["CVN"],Ind["CVG"],
             Ind["HP"],Ind["C"])
    ind <- c(ind,TL=TT$TL[specs])
    ind <- c(ind,OI=TT$OI[specs])
    
    # and add them to the results matrices
    Indices <- rbind(Indices,Ind)
    indices <- rbind(indices,ind)
  }
  list(Full=Indices,   # all calculated indices
       Sub=indices)    # selected set
}
# import LIM results for each food web

S1.lim$NComponents

S1.ENA<-Indices(S1.lim, S1.X, ext=c("gpp","res","los"),dead=c(7,8), specs=c(5,6))



##calculation of detritivory/herbivory indice##
  
  DH <-numeric(nrow(S1.X))

  S1.X<- res_S1_300000
for (j in 1:nrow(S1.X)){
  DH[j]<- (sum(S1.X[j,33])/sum(S1.X[j,c(7,8,13,14,18)]))
}
DH
DH<-sample(DH,size=300000,replace=FALSE)
DH1<-as.data.frame(DH)
summary(DH1)
saveRDS(DH1,"C:/Users/DELL/.../DH_S11.rds")

#the same calculation of ENA and DH was repeated for the other stations, each time replacing S1.X  and S1.lim by S2.X / S2.lim then S3.X/S3/lim then S4.X/S4.lim
  ```

### Plot ENA_indices
```R
ENA.S1<-gather(ENA.S1, indices)%>%
  mutate( web="S1")
DH1<-gather(DH1, indices)%>%
  mutate( web="S1")

ENA.S2<-gather(ENA.S2, indices)%>%
  mutate( web="S2")
DH2<-gather(DH2, indices)%>%
  mutate( web="S2")
#ENAS2<- rbind(ENA.S2,DH2)

ENA.S3<-gather(ENA.S3, indices)%>%
  mutate( web="S3")
DH3<-gather(DH3, indices)%>%
  mutate( web="S3")
#ENAS3<- rbind(ENA.S3,DH3)

ENA.S4<-gather(ENA.S4, indices)%>%
  mutate( web="S4")
DH4<-gather(DH4, indices)%>%
  mutate( web="S4")
#ENAS4<- rbind(ENA.S4,DH4)


IND_ENA<-bind_rows(ENA.S1,ENA.S2,ENA.S3,ENA.S4)
IND_ENA2<-bind_rows(DH1,DH2,DH3,DH4)
IND_ENA22<-bind_rows(IND_ENA,IND_ENA2)



X<-filter(IND_ENA22_bac,indices %in% c("TST", "AMI", "APL", "FCI","ACratio","DH"))%>%
  mutate(indices=as.factor(indices))

p<-ggplot(data=X,aes(x=web,y=value))
p+geom_boxplot()+theme_classic()+facet_wrap(~indices,scale="free")+
  theme(strip.text.x=element_text(size=16), axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), axis.title =element_text(size=16))+
  labs(x="Web",y="values")
```
![ENA_indices](https://github.com/Chkili/article_indicators/blob/main/ENA_indices.jpg)
 
 Spatial variation of ENA indices calculated for the planktonic food webs in four stations of the Gulf of Gabès. Total system throughput (TST; mg C m−2 d−1) (A), relative ascendancy (A/C; %) (B), average mutual information (AMI; bits) (C), average path length (APL) (D), cycling index (FCI; %) (E), and detritivory to herbivory (D/H) (F).]

## Calculation of Delta Cliff

The Cliff’s δ test was used to statistically test indices differences between models.

```R
S1<-readRDS("C:/Users/..../S1_RES.ENA.rds")
S2<-readRDS("C:/Users/..../S2_RES.ENA.rds")
S2<-readRDS("C:/Users/..../S2_RES.ENA.rds")
S1<-readRDS("C:/Users/..../S1_RES.ENA.rds")

DH1<-readRDS("C:/Users/..../DH1.rds")
DH2<-readRDS("C:/Users/..../DH2.rds")
DH3<-readRDS("C:/Users/..../DH3.rds")
DH4<-readRDS("C:/Users/..../DH4.rds")

df_S1 <- S1$Sub[,c(1,6,7,15,16)]
df_S1 <- as.data.frame(df_S1)

colnames(S2$Sub)

df_S2 <- S2$Sub[,c(1,6,7,15,16)]
df_S2 <- as.data.frame(df_S2)

r2=NULL

ENAnames <- colnames(df_S2)

for (c in 1 : 5){ # choisir les colonnes que tu veux (les ENA)
  
  x <- df_S1[1:100,c] # mod�le de ref du cliff, tu compare donc S2 � S2
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
  assign(paste0("S1_S2_ENA_Cliff"), r2) # la tu cr�er une dataframe avec tes resultats
}

saveRDS(S1_S2_ENA_Cliff, file="C:/Users/..../S1_S2_ENA_Cliff.rds")
```
## Calculation of Multiple factor analysis (MFA)

A multiple factor analysis (MFA) was performed to identify the interrelationships between different ecological indicators (food web typology ratios and ENA indices) as well as environmental variables (inorganic and organic nutrients) and some of the calculated fluxes (GPP of PIC, NAN and MIC, bacterial production and sinking of NAN, MIC and MET). 

#### Environmental data
[Environmental data](https://github.com/Chkili/article_indicators/blob/main/envi_golfe.rds)
<br>

``` R
calculation of 3000 randomly selected values LIM solutions 

filesr   = list.files("C:/Users/....", pattern = "_variables.rds", full.names = TRUE)


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
```
![MFA](https://github.com/Chkili/article_indicators/blob/main/MFA..jpg)
Multiple factor analysis (MFA) ordination diagram showing the relationships between ecological indicators (food web typology ratios and ENA indices), environmental variables (inorganic and organic nutrients: Ninorg, Pinorg, Siinorg, Norg, Porg) and carbon flows [GPP of PIC (GPP->PIC), NAN (GPP->NAN) and MIC (GPP->MIC), bacterial production (DOC->BAC) and sinking of NAN (NAN->LOS), MIC (MIC->LOS) and MET (MET->LOS)]. 

  
  
 



  
 




