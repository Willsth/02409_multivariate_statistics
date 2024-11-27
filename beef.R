Beef <- data.frame("pH" = c(1,rep(0,17)),
                "Water" = c(0.09,1,rep(0,16)),
                "Protein"  = c(0.28,-0.4,1,rep(0,15)),
                "EtherExt" = c(-0.28,-0.16,-0.56,1,rep(0,14)),
                "Hydroxy" = c(-0.33,-0.08,-0.55,0.59,1,rep(0,13)),
                "CollaSol" = c(-0.08,-0.01,-0.03,0.05,0.16,1,rep(0,12)),
                "Lightn" = c(-0.02,0.03,0.34,-0.31,-0.48,-0.02,1,rep(0,11)),
                "Hue" = c(-0.33,-0.23,-0.47,0.4,0.62,-0.03,-0.21,1,rep(0,10)),
                "DripLoss" = c(0.01,0.18,-0.07,0.08,-0.12,-0.1,0.25,-0.13,1,rep(0,9)),
                "CookLoss" = c(-0.38,0.15,-0.64,0.44,0.66,-0.01,-0.45,0.65,0.03,1,rep(0,8)),
                "WBshear" = c(-0.26,-0.01,-0.63,0.42,0.72,-0.03,-0.55,0.67,-0.11,0.73,1,rep(0,7)),
                "Appear" = c(0.1,-0.003,0.25,-0.42,-0.33,-0.19,0.35,0.07,0.02,-0.18,-0.28,1,rep(0,6)),
                "EaseSink" = c(0.17,-0.16,0.27,-0.11,-0.26,0.01,0.19,-0.19,-0.02,-0.36,-0.38,0.24,1,rep(0,5)),
                "Friabil" = c(0.1,-0.17,0.2,-0.09,-0.22,0.06,0.19,-0.1,-0.03,-0.31,-0.32,0.27,0.93,1,rep(0,4)),
                "Residue"  = c(0.08,-0.19,0.23,-0.13,-0.24,-0.02,0.2,-0.1,-0.02,-0.32,-0.33,0.33,0.91,0.94,1,rep(0,3)),
                "InJuice" = c(0.08,-0.08,0.03,-0.004,-0.05,0.05,-0.06,-0.03,-0.13,-0.12,-0.1,0.16,0.69,0.72,0.72,1,rep(0,2)),
                "SusJuice" = c(0.01,-0.09,-0.004,-0.01,-0.02,0.03,-0.02,0.08,-0.15,-0.07,-0.03,0.24,0.66,0.7,0.7,0.93,1,rep(0,1)),
                "OvAcc" = c(0.13,-0.13,0.21,-0.09,-0.22,0.07,0.22,-0.14,-0.01,-0.34,-0.37,0.31,0.92,0.92,0.91,0.8,0.79,1)
)

Beef = Beef + t(Beef) - diag(1,18,18)

# Recommended Packages:
library(CCP)
install.packages("geigen")  # Install the package if not already installed
library(geigen)

Exx = as.matrix(Beef[1:11,1:11])
Eyx = as.matrix(Beef[12:18,1:11])
Exy = as.matrix(Beef[1:11,12:18])
Eyy = as.matrix(Beef[12:18,12:18])
invExx = solve(Exx)
invEyy = solve(Eyy)

#Calculating the Canonical correlations:
Cancorr = geigen(Eyx%*%invExx%*%Exy,Eyy,symmetric = TRUE)
values = sort(Cancorr$values,decreasing = TRUE)

# E is the residual variation after having predicted Y by means of X
H = Eyx%*%invExx%*%Exy
E = Eyy - Eyx%*%invExx%*%Exy
invE = solve(E)
Ev <- eigen(invE%*%H)
var = Ev$values

# Eigenvalues, Proportion and Cumulative proportion of Variance:
varPC <- var/sum(var)
cumu = c(1:7)
for (i in 1:7){
  cumu[i] = sum(varPC[1:i])
}
results <- data.frame("CanCor" = sqrt(values),"Squared CanCor" = values,"eigenvalues"=var,"proportion"=varPC,
                      "cumulative" = cumu)
print("Table with information about the Canonical Correlations:")
results