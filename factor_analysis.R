Heptathlon = data.frame("Name" = c("Joyner_Kersee,(USA)","John,(GDR)","Behmer,(GDR)","Choubenkova,(URS)","Sablovskaite,(URS)","Schulz,(GDR)","Fleming,(AUS)","Greiner,(USA)","Lajbnerova,(CZE)","Bouraga,(URS)","Wijnsma,(HOL)","Dimitrova,(BUL)","Scheider,(SWI)","Braun,(FRG)","Ruotsalainen,(FIN)","Yuping,(CHN)","Hagger,(GB)","Brown,(USA)","Mulliner,(GB)","Hautenauve,(BEL)","Kytola,(FIN)","Geremias,(BRA)","Hui-Ing,(TAI)","Jeong-Mi,(KOR)","Launa,(PNG)"),
                      "Points"= c(7291,6897,6858,6540,6456,6411,6351,6297,6252,6232,6205,6171,6157,6109,6101,6087,5975,5972,5746,5734,5686,5508,5290,5289,4566),
                      "Hurdles" = c(12.69,12.85,13.20,13.51,13.61,13.75,13.38,13.55,13.63,13.25,13.75,13.24,13.85,13.71,13.79,13.93,13.47,14.07,14.39,14.04,14.31,14.23,14.85,14.53,16.42),
                      "High Jump"=c(1.86,1.80,1.83,1.74,1.80,1.83,1.80,1.80,1.83,1.77,1.86,1.80,1.86,1.83,1.80,1.86,1.80,1.83,1.71,1.77,1.77,1.71,1.68,1.71,1.50),
                       "Shot" = c(15.80,16.23,14.20,14.76,15.23,13.50,12.88,14.13,14.28,12.62,13.01,12.02,11.58,13.16,12.32,14.21,12.75,12.69,12.68,11.81,11.66,12.95,10.00,10.83,11.78),
                       "Run200" = c(22.56,23.65,23.10,23.93,23.92,24.65,23.59,24.48,24.86,23.59,25.03,23.59,24.87,24.78,24.61,25.00,25.47,24.83,24.92,25.61,25.69,25.50,25.23,26.61,26.16),
                       "Longjump" = c(7.27,6.71,6.68,6.32,6.25,6.33,6.37,6.47,6.11,6.28,6.34,6.19,6.05,6.12,6.08,6.40,6.34,6.13,6.10,5.99,5.75,5.50,5.47,5.50,4.88),
                       "Javelin" = c(45.66,42.56,44.54,47.46,42.78,42.82,40.28,38.00,43.30,39.06,37.86,37.62,47.50,44.58,45.44,38.60,35.78,44.34,37.76,35.68,39.48,39.64,39.14,39.26,46.38),
                       "Run800" = c(128.51,126.14,124.20,127.90,132.24,125.79,132.54,133.65,136.05,134.74,131.49,135.73,134.93,142.82,137.06,146.67,138.48,146.43,138.02,133.90,133.35,144.02,137.30,139.17,163.43))



Hep <-Heptathlon[,-c(1,2)]
Eigs <- eigen(cor(Hep))

var <- Eigs$values
varProp <- var/sum(var)
cumu = rep(NA,length(var))
for (i in 1:length(var)){
  cumu[i] = sum(varProp[1:i])
}
results <- data.frame("eigenvalues"=var,"proportion"=varProp,
                      "cumulative" = cumu)
print("Eigenvalues of the Correlation Matrix:")
results     


# load principal package

library(psych)
#### Factor analysis 2 Factors ####
#unrotated factors
fa <- principal(cor(Hep),nfactors = 2,rotate = "none")
fal <- fa$loadings[,1:2]
facom <- fa$communality

#rotated factors:
rfa <-principal(cor(Hep),nfactors = 2,rotate = "varimax",scores = T)
rfal <- rfa$loadings[,1:2]
rfacom <- rfa$communality

# The 1st Factor is correct as it comes directly from pca's Evals and Evecs, but rotation...?
Factors <- data.frame("FA" = fal,"Rot_FA" = rfal)
print("Non rotated and rotated factors:")
Factors

#Communality:
Commu <- data.frame("comun" = facom, "rot_comun" = rfacom)
print("Final Communality Estimates (non rotated and rotated):")
Commu

#Variance explained:
var <- c(sum(fal[,1]^2),sum(fal[,2]^2))
# var can also be found with fa$values[1:2], but not available for rotated factors
rvar <- c(sum(rfal[,1]^2),sum(rfal[,2]^2))
Var <- data.frame("var" = var,"rot_var" = rvar)
print("Variance Explained by Each Factor:")
Var

sum(1/7*Var[,1])*100
sum(1/7*Var[,2])*100

Factors$FA.PC1[1]^2+Factors$FA.PC2[1]^2
Factors$Rot_FA.RC1[1]^2+Factors$Rot_FA.RC2[1]^2


par(mfrow=c(1,2))
circle = seq(-3.2,3.2,by=0.1) 

#Plot for Factors 1&2
plot(0,0,xlim = c(-1.2,1.2),ylim = c(-1.2,1.2),xlab = "Factor 1",
     ylab = "Factor 2",main = "Initial Factor Pattern")
points(1*cos(circle),1*sin(circle),type='l')
arrows(c(rep(0,7)),c(rep(0,7)),fal[,1],fal[,2],length = 0.1)
text(fal[,1],fal[,2]+0.1,names(Hep),cex = 0.7)
grid()

#Plot for rotated Factors 1&2
plot(0,0,xlim = c(-1.2,1.2),ylim = c(-1.2,1.2),xlab = "Rotated Factor 1",
     ylab = "Rotated Factor 2",main = "Rotated Factor Pattern")
points(1*cos(circle),1*sin(circle),type='l')
arrows(c(rep(0,7)),c(rep(0,7)),rfal[,1],rfal[,2],length = 0.1)
text(rfal[,1],rfal[,2]+0.1,names(Hep),cex = 0.7)
grid()


fa3 <- principal(cor(Hep),nfactors = 3,rotate = "none")
fa3l <- fa3$loadings[,1:3]
fa3com <- fa3$communality

#rotated:
rfa3 <- principal(cor(Hep),nfactors = 3,rotate = "varimax") 
rfa3l <- rfa3$loadings[,1:3]
rfa3com <- rfa3$communality

#Factors
Factors3 <- data.frame("FA" = fa3l,"Rot FA" = rfa3l)
print("Non rotated and rotated factors:")
Factors3

#Communality:
Commu3 <- data.frame("comun" = fa3com, "rot comun" = rfa3com)
print("Final Communality Estimates (non rotated and rotated):")
Commu3

#Variance explained:
rvar3 = c(sum(rfa3l[,1]^2),sum(rfa3l[,2]^2),sum(rfa3l[,3]^2))
var3 <- c(sum(fa3l[,1]^2),sum(fa3l[,2]^2),sum(fa3l[,3]^2))
Var3 <- data.frame("var" = var3[1:3],"rotvar" = rvar3)
print("Variance Explained by Each Factor:")
Var3

# Plots for factor analysis with 3 factors:
par(mfrow = c(1,2))
circle = seq(-3.2,3.2,by=0.1)

# Different combinations of plots
ij = matrix(c(1,1,2,2,3,3),ncol=2)

Names = c("Hdls","HgJmp","Shot","R200","LgJmp","Jvln","R800")
for (i in 1:3){
  l = ij[i,1]
  k = ij[i,2]
  
  #Plot for the Factors
  plot(0,0,xlim = c(-1.2,1.2),ylim = c(-1.2,1.2),xlab = paste0("Factor ",l),
       ylab = paste0("Factor ",k),main = "Initial Factor Pattern")
  points(1*cos(circle),1*sin(circle),type='l')
  arrows(c(rep(0,7)),c(rep(0,7)),fa3l[,l],fa3l[,k],length = 0.1)
  text(fa3l[,l],fa3l[,k]+0.1,Names,cex = 0.7)
  grid()
  
  #Plot for rotated Factors
  plot(0,0,xlim = c(-1.2,1.2),ylim = c(-1.2,1.2),xlab = paste0("Factor ",l),
       ylab = paste0("Factor ",k),main = "Rotated Factor Pattern")
  points(1*cos(circle),1*sin(circle),type='l')
  arrows(c(rep(0,7)),c(rep(0,7)),rfa3l[,l],rfa3l[,k],length = 0.1)
  text(rfa3l[,l],rfa3l[,k]+0.1,Names,cex = 0.7)
  grid()
}