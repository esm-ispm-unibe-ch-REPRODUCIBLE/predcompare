expit<-function(x){ exp(x)/(1 + exp(x))}

simulate.data.cont<-function(Npat){
  ### simulate covariates
  x1x2<-mvrnorm(n = Npat, c(0,0), Sigma=matrix(c(1,0.2,0.2,1),2,2))
  simdat <- data.frame(x1=x1x2[,1]) # continuous covariate 1
  simdat$x2<-x1x2[,2] # continuous covariate 2
  x3.lin <-1+simdat$x1-0.2* simdat$x2+rnorm(Npat, 0, 1)
  simdat$x3<-rbinom(Npat, 1, prob = expit(x3.lin) )# binary covariate 13
  
  simdat$x4<-rbinom(Npat, 1, prob = 0.2 )
  simdat$x5<-rnorm(Npat, 0, 1 )
  ### simulate treatment assignment
  pt <- 0.5 ## randomized
  simdat$t <- rbinom(Npat, 1, prob = pt)
  ### outcome at control treatment
  simdat$y.control<- with(simdat, 1+0.3*x1+0.5*x2+0.2*x3+0.3*x4+0.4*x5+
                            0.2*x1^2+0.1*log(abs(2*x1+1))+0.2*x1*x2+0.2*x1*x2*x3+ 
                            rnorm(Npat, 0, 0.2))
  ### treatment benefit
  simdat$benefit<-with(simdat, -0.1+0.2*x1+0.1*x2+0.1*x3+
                         0.1*x1*x2+0.1*x1*x3+0.1*sqrt(x1^2+x2^2)+
                         rnorm(Npat, 0, 0.3))
  simdat$y.active<-with(simdat,y.control+benefit)
  simdat$y.observed<-with(simdat,y.control+t*benefit)
  return(list(dat=simdat))}


