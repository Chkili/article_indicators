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
