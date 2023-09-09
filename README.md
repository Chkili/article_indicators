# Journal of Marine Systems 
This R script comes with the following research article:
It includes (i) a procedure for calculation of the unknown fluxes using the LIM-MCMC method based on the mirror technique defined by Meersche et al. (2009) estimates each unknown flux, (ii) a procedure for estimating ENA indicators and typology ratios and (iv) procedure to create a Multiple Factor Analysis (MFA) Ordination Diagram showing the relationships between ecological indicators, environmental variables and carbon fluxes.

## Data_matrix
[Station1-S1](https://www.google.com)
[Station2-S2](https://www.google.com)
[Station3-S3](https://www.google.com)
[Station4-S4](https://www.google.com)

## Load required libraries

`library(LIM)
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
library(stringr)`

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



