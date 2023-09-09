setwd("C:/Users/DELL/...") #to set which directory you want to work in: source path
getwd() #to check the right path
S1.lim <- Setup("Golfe_S1.txt") #download declaration file and transformation into ABGH matrix
pars <- Ldei(S1.lim) #pars it's to calculate the deterministic (least squares) (parsimonious) solution (a single value for a flow instead of a set of solutions)
webranges<- Xranges(S1.lim) #xranges it's to finds the possible ranges ([min,max]) for each unknown.
S1.X<-Xsample(S1.lim,iter=300000,jmp=10) # calculation of random solution jump=10 iter=300000
Summary(S1.X)