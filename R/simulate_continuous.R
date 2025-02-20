#' Simulate data for a prediction model of a continuous outcome
#'
#' This function generates a dataframe with 6 patient covariates and a continuous outcome 
#' simulated from a model that uses the covariates. 
#' @param Npat Number of patients to simulate.
#' @return The function returns a dataframe with: 
#' 
#' x1, x2, x3, x4, x5, x6= patient covariates.
#' 
#' t= treatment assignment (0 for control, 1 for active).
#' 
#' y.control= the outcome if the patient takes the control treatment. 
#' 
#' y.active= the outcome if the patient takes the active treatment. 
#' 
#' benefit= the treatment benefit, i.e. y.active-y.control.
#' 
#' y.observed= the observed outcome. 
#' @examples 
#' dat1=simcont(100)$dat
#' head(dat1)
#' @export
simcont=function(Npat){
  library(MASS)
  ### simulate covariates
  expit=function(x){ exp(x)/(1 + exp(x))}
  # observed true predictors
  x1x2<-mvrnorm(n = Npat, c(0,0), Sigma=matrix(c(1,0.2,0.2,1),2,2))
  simdat <- data.frame(x1=x1x2[,1])
  simdat$x2<-x1x2[,2] 
  simdat$x3<-rbinom(Npat, 1, prob = 0.2)
  simdat$x4<-rnorm(Npat, 0, 1 )
  ##  nuisance parameters
  simdat$x5<-rnorm(Npat, 0, 1 )
  simdat$x6<-rnorm(Npat, 0, 1 )
  
  ### simulate treatment assignment
  pt <- 0.5 ## randomized
  simdat$t <- rbinom(Npat, 1, prob = pt)
  ### outcome at control treatment
  simdat$y.control<- with(simdat, 1+0.3*x1+0.5*x2+0.2*x3+0.3*x4+
                            0.2*x1^2+0.2*log(abs(2*x1+1))+0.2*x1*x2+0.2*x1*x2*x3+ 
                            rnorm(Npat, 0, 0.2))
  ### treatment benefit
  simdat$benefit<-with(simdat, -0.1+0.2*x1+0.1*x2+0.1*x3+0.2*x4+
                         0.1*x1*x2+0.1*x1*x3+0.1*sqrt(x1^2+x2^2)+
                         rnorm(Npat, 0, 0.2))
  simdat$y.active<-with(simdat,y.control+benefit)
  simdat$y.observed<-with(simdat,y.control+t*benefit)
  return(list(dat=simdat))}


