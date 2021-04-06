rm(list = ls()) # Removing R's memory

# Where you want to save the figures (the path)

whereAmI = "C:/Users/jtomal/Desktop/Metron/Metron/"

# Reading data from external source

cdata = read.csv("C:/Users/jtomal/Desktop/Metron/Metron/MATH_1070.csv" , header=TRUE)

# Checking the first few rows of the data

head(cdata)

# Checking the dimension of the data (number of rows and columns)

dim(cdata)

# The last row contains the covid19 (1 = After March 15, 0 = Before March 15) variable

covid19 <- c(cdata[dim(cdata)[[1]],])
covid19

# The data now only contains marks

cdata <- cdata[-dim(cdata)[[1]],]
head(cdata)

# The number of students in the course

m <- dim(cdata)[[1]]
m

# The number of evaluation components (dimension)

d <- dim(cdata)[[2]]
d

# Checking if any student disappeared either after or before the COVID effect

disappeared.student <- NULL

for(i in 1:m){
  if(all(is.na(cdata[i, (covid19 == 1)])) || all(is.na(cdata[i, (covid19 == 0)]))){disappeared.student <- c(disappeared.student, i)}
}

disappeared.student

# Proportion of disappeared students

PD <- (length(disappeared.student)/m)*100
PD

# Removing disappeared students' rows from the data

cdata <- cdata[-disappeared.student,]
dim(cdata)

# Number of remaining students for which we have marks for analysis

m <- dim(cdata)[[1]]
m

# Number of evaluation components (dimension)

d <- dim(cdata)[[2]]
d

# which marks are observed (O = 1) and which are missing (O = 0)

Y.miss <- cdata     # This data contains both observed and missing marks
Y.miss
O <- 1*(!is.na(Y.miss))           # Index for missing: O = 1 for Observed, O = 1 for missing
O

((length(O)-sum(O))/length(O))*100  # Proportion of missing marks


# Which rows contain missing data? The object 'miss' gives the row numbers

miss <- NULL
for(i in 1:m){
  for(j in 1:d){
    if(is.na(Y.miss[i,j])){miss <- c(miss, i)}
  }
}
miss <- sort(unique(miss))
miss

# Y.full is complete data where the missing values are replaced by the average of available data

Y.full <- Y.miss
for(j in 1:d)
{
  Y.full[is.na(Y.full[,j]), j]<-mean(Y.full[,j], na.rm=TRUE)
}
Y.full
row.names(Y.full) <- colnames(Y.full) <- NULL

# Create the data in list format

one.vec <- rep(1, d)
design.mat <- cbind(one.vec, covid19)
row.names(design.mat) <- NULL
colnames(design.mat) <- c("x.0", "x.covid19")
design.mat

Y<-list() ; X<-list()
for(j in 1:m) 
{
  Y[[j]]<-unlist(Y.full[j,])
  X[[j]]<- matrix(unlist(cbind( rep(1, d), covid19  )), ncol = 2, byrow = F)
}

#### OLS (Ordinary list-squares) fits: Linear regression fit to each student's marks

S2.LS <- BETA.LS <- NULL
for(i in 1:m) {
  fit <- lm(Y[[i]]~-1+X[[i]] )
  BETA.LS <- rbind(BETA.LS, c(fit$coef)) 
  S2.LS <- c(S2.LS, summary(fit)$sigma^2) 
}



#### Hierarchical regression model

## mvnormal simulation
rmvnorm<-function(n,mu,Sigma)
{ 
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}

## Wishart simulation
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}


## Setup

p <- dim(X[[1]])[2]
theta <- mu0 <- apply(BETA.LS,2,mean)
nu0 <- 1 ; s2 <- S2.LS
s20<- mean(S2.LS)
nu0s <- 1:5000
alpha1 <-1 ; alpha2 <- 1/100 ; alpha3 <- 1

eta0 <- p+2 ; Sigma <- S0<- L0<- cov(BETA.LS) ; BETA <- BETA.LS
THETA.b <- S2.b <- S20.b <- NU0 <- NULL
iL0 <- solve(L0) ; iSigma <- solve(Sigma)
Sigma.ps <- matrix(0,p,p)
SIGMA.PS <- NULL
BETA.ps.a <- BETA.ps.b <- NULL #BETA*0
mu0[2]+c(-1.96,1.96)*sqrt(L0[2,2])
B <- 1000 # Burn-in period
R <- 11000

set.seed(1)

## MCMC (Markov Chain Monte Carlo)
for(r in 1:R) {

  #sample new s20 (the overall variance)
  s20<-rgamma(1, shape = alpha1+(m*nu0)/2, rate = alpha2+(nu0/2)*sum(1/s2))
  
  # sample nu0 (the strength of the overall variance)
  
  lpnu0<- .5*nu0s*m*log(s20*nu0s/2)-m*lgamma(nu0s/2)+(nu0s/2-1)*sum(log(1/s2)) -
    nu0s*s20*sum(1/s2)/2   - alpha3*nu0s
  nu0<-sample(nu0s,1,prob=exp( lpnu0-max(lpnu0)) )
  

  ##update s2
  for(j in 1:m) {
    RSS <- 0
    RSS <- (Y[[j]]-X[[j]]%*%BETA[j,] )^2
    s2[j] <-1/rgamma(1, shape = (nu0+ d)/2, rate = (nu0*s20+RSS)/2 )
  }
  ##
  
    ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma + t(X[[j]])%*%X[[j]]/s2[j] )
    Ej<-Vj%*%( iSigma%*%theta + t(X[[j]])%*%Y[[j]]/s2[j] )
    BETA[j,]<-rmvnorm(1,Ej,Vj) 
  } 
  ##
  
  ##update theta
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum))
  theta<-t(rmvnorm(1,mum,Lm))
  ##
  
  ##update Sigma
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA-mtheta)%*%(BETA-mtheta) ) )
  ##
  

  ###update missing data
  for(j in miss)
  { ESigma <- matrix(0, ncol = d, nrow = d)
    diag(ESigma) <- s2[j]
    XBETA <- X[[j]]%*%BETA[j,]
    YSigma <- (X[[j]] %*% solve(iSigma) %*% t(X[[j]])) + ESigma
    b <- ( O[j,]==0 )
    a <- ( O[j,]==1 )
    iSa<- solve(YSigma[a,a])
    beta.j <- YSigma[b,a]%*%iSa
    s2.j   <- YSigma[b,b] - YSigma[b,a]%*%iSa%*%YSigma[a,b]
    theta.j<- XBETA[b] + beta.j%*%(t(Y.miss[j,a])-XBETA[a])
    Y.miss[j,b] <- rmvnorm(1,theta.j,s2.j )
  }
  
  Y<-list() ; X<-list()
  for(j in 1:m) 
  {
    Y[[j]]<-unlist(Y.miss[j,])
    X[[j]]<- matrix(unlist(cbind( rep(1, d), covid19  )), ncol = 2, byrow = F)
  }
  
  ##store results
  if((r%%10==0) & (r > B)) 
  { 
    cat("Iteration #: ", r, "; Overall error Variance:", s20,"\n")
    NU0 <- c(NU0, nu0)
    S20.b <- c(S20.b, s20)
    S2.b<-rbind(S2.b,s2)
    THETA.b<-rbind(THETA.b,t(theta))
    SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
    BETA.ps.a <- rbind(BETA.ps.a, c(BETA[,1])) 
    BETA.ps.b <- rbind(BETA.ps.b, c(BETA[,2])) 
    }
  ##
  
}


# Plot the CIs for \theta (overall marks in a course)

par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))
df <- data.frame(x =1:2,
                 F = c(quantile(THETA.b[,1], probs = 0.50), quantile(apply(THETA.b, 1, sum), probs = 0.50)),
                 L = c(quantile(THETA.b[,1], probs = 0.025), quantile(apply(THETA.b, 1, sum), probs = 0.025)),
                 U = c(quantile(THETA.b[,1], probs = 0.975), quantile(apply(THETA.b, 1, sum), probs = 0.975)))

require(ggplot2)
ggplot(df, aes(x = x, y = F)) +
  labs(y = "Marks (%)", x = "COVID19 Effects") +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Before","After")) +
  theme(axis.text.x = element_text(color = "black"),
  axis.ticks.x = element_line(color = c("black", "black"),
  size = c(1, 1)))
#dev.print(pdf,paste(whereAmI,"MATH_1070_Theta_CI.pdf",sep=""),width=8,height=8)


# Create latex table for overall marks in a course

TB <- cbind(c(quantile(THETA.b[,1], probs = c(0.025, 0.50, 0.975)), sd(THETA.b[,1])),
        c(quantile(apply(THETA.b, 1, sum), probs = c(0.025, 0.50, 0.975)), sd(apply(THETA.b, 1, sum))))
row.names(TB) <- c("Lower Limit", "Median", "Upper Limit", "SD")
colnames(TB) <- c("Before COVID", "After COVID")
TB
#install.packages("xtable")
library(xtable)
xtable(TB, digits = 3)


# Plot the distribution/density for \theta

plot(density(apply(THETA.b, 1, sum),adj=2),xlim=range(c(THETA.b[,1], apply(THETA.b, 1, sum))), ylim = c(0, 0.58), 
     main="",xlab="",ylab="",  lwd=4, cex.lab=1.7, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
lines(density(THETA.b[,1],adj=2),col="gray",lwd=2)
legend( 56, 0.58 ,legend=c( expression(theta[0]),expression(theta[0] + theta[1])), 
        lwd=c(2,2),col=c("gray", "black"),cex=1.7,bty="n") 
#dev.print(pdf,paste(whereAmI,"MATH_1070_Theta_Dist.pdf",sep=""),width=8,height=8)



# Plot the distribution/density of \beta

i <- 5
plot(density(BETA.ps.a[,i] + BETA.ps.b[,i],adj=2),xlim=range(c(BETA.ps.a[,i], BETA.ps.a[,i] + BETA.ps.b[,i])), ylim=c(0, 0.58), 
     main="",xlab="",ylab="",  lwd=4, cex.lab=1.7, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
lines(density(BETA.ps.a[,i],adj=2),col="gray",lwd=2)
legend( 60, 0.58 ,legend=c( expression(beta[0]),expression(beta[0] + beta[1])), 
        lwd=c(2,2),col=c("gray", "black"),cex=1.7,bty="n") 
#dev.print(pdf,paste(whereAmI,paste(paste("Math_1070_Dist_Student", i, sep = "_"), "pdf", sep = "."),sep=""),width=8,height=8)

# Create latex table for ith student

TB <- cbind(c(quantile(BETA.ps.a[,i], probs = c(0.025, 0.50, 0.975)), sd(BETA.ps.a[,i])),
            c(quantile(BETA.ps.a[,i] + BETA.ps.b[,i], probs = c(0.025, 0.50, 0.975)), sd(BETA.ps.a[,i] + BETA.ps.b[,i])))
row.names(TB) <- c("Lower Limit", "Median", "Upper Limit", "SD")
colnames(TB) <- c("Before COVID", "After COVID")
TB
#install.packages("xtable")
library(xtable)
xtable(TB, digits = 3)


i <- 14
plot(density(BETA.ps.a[,i] + BETA.ps.b[,i],adj=2),xlim=range(c(BETA.ps.a[,i], BETA.ps.a[,i] + BETA.ps.b[,i])), ylim=c(0, 0.58), 
     main="",xlab="",ylab="",  lwd=4, cex.lab=1.7, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
lines(density(BETA.ps.a[,i],adj=2),col="gray",lwd=2)
legend( 75, 0.58 ,legend=c( expression(beta[0]),expression(beta[0] + beta[1])), 
        lwd=c(2,2),col=c("gray", "black"),cex=1.7,bty="n") 
#dev.print(pdf,paste(whereAmI,paste(paste("Math_1070_Dist_Student", i, sep = "_"), "pdf", sep = "."),sep=""),width=8,height=8)

# Create latex table for ith student

TB <- cbind(c(quantile(BETA.ps.a[,i], probs = c(0.025, 0.50, 0.975)), sd(BETA.ps.a[,i])),
            c(quantile(BETA.ps.a[,i] + BETA.ps.b[,i], probs = c(0.025, 0.50, 0.975)), sd(BETA.ps.a[,i] + BETA.ps.b[,i])))
row.names(TB) <- c("Lower Limit", "Median", "Upper Limit", "SD")
colnames(TB) <- c("Before COVID", "After COVID")
TB
#install.packages("xtable")
library(xtable)
xtable(TB, digits = 3)

# Figure for individual student

i <- 27
plot(density(BETA.ps.a[,i] + BETA.ps.b[,i],adj=2),xlim=range(c(BETA.ps.a[,i], BETA.ps.a[,i] + BETA.ps.b[,i])), ylim=c(0, 0.58), 
     main="",xlab="",ylab="",  lwd=4, cex.lab=1.7, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
lines(density(BETA.ps.a[,i],adj=2),col="gray",lwd=2)
legend( 77, 0.58 ,legend=c( expression(beta[0]),expression(beta[0] + beta[1])), 
        lwd=c(2,2),col=c("gray", "black"),cex=1.7,bty="n") 
#dev.print(pdf,paste(whereAmI,paste(paste("Math_1070_Dist_Student", i, sep = "_"), "pdf", sep = "."),sep=""),width=8,height=8)


TB <- cbind(c(quantile(BETA.ps.a[,i], probs = c(0.025, 0.50, 0.975)), sd(BETA.ps.a[,i])),
            c(quantile(BETA.ps.a[,i] + BETA.ps.b[,i], probs = c(0.025, 0.50, 0.975)), sd(BETA.ps.a[,i] + BETA.ps.b[,i])))
row.names(TB) <- c("Lower Limit", "Median", "Upper Limit", "SD")
colnames(TB) <- c("Before COVID", "After COVID")
TB
#install.packages("xtable")
library(xtable)
xtable(TB, digits = 3)