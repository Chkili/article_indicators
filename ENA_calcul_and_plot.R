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

saveRDS(S1.ENA, file="C:/Users/DELL/..../S1_ENA.rds")




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

#the same calculation of ENA and DH was repeated for the other stations, 
#each time replacing S1.X  and S1.lim by S2.X / S2.lim then S3.X/S3/lim then S4.X/S4.lim

![](https://github.com/Chkili/article_indicators/blob/main/ENA_indices.jpg)
