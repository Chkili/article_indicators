# Journal of Marine Systems 
This R script comes with the following research article:
It includes (i) a procedure for calculation of the unknown fluxes using the LIM-MCMC method based on the mirror technique defined by Meersche et al. (2009) estimates each unknown flux, (ii) a procedure for estimating ENA indicators and typology ratios and (iv) procedure to create a Multiple Factor Analysis (MFA) Ordination Diagram showing the relationships between ecological indicators, environmental variables and carbon fluxes.


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