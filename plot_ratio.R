ratio.S1<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/Rratio/ratioo/ratioS1/S1_ratio.rds")
ratio.S2<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/Rratio/ratioo/ratioS2/S2_ratio.rds")
ratio.S3<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/Rratio/ratioo/ratioS3/S3_ratio.rds")
ratio.S4<-readRDS("C:/Users/DELL/Documents/Ancien/modelisation gabes/Rratio/ratioo/ratioS4/S4_ratio.rds")

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

