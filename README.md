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
library("FactoMineR")
library("factoextra")
library(qgraph)
library(igraph)
library("ggplot2")
library(effsize)
library(stringr)
```

## Calculation of the unknown fluxes
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

![ENA_indices](https://github.com/Chkili/article_indicators/blob/main/ENA_indices.jpg)
 
 Spatial variation of ENA indices calculated for the planktonic food webs in four stations of the Gulf of Gabès. Total system throughput (TST; mg C m−2 d−1) (A), relative ascendancy (A/C; %) (B), average mutual information (AMI; bits) (C), average path length (APL) (D), cycling index (FCI; %) (E), and detritivory to herbivory (D/H) (F).]
  
 




