#' Simulate data for a prediction model of a binary outcome
#'
#' This function generates a dataframe with 6 patient covariates and a binary outcome simulated from
#' a model that uses the covariates.  
#' @param Npat Number of patients to simulate.
#' @return The function returns a dataframe with: 
#' 
#' x1, x2, x3, x4, x5, x6= patient covariates.
#' 
#' t= treatment assignment (0 for control, 1 for active).
#' 
#' logit.control= the logit of the probability of an outcome in the control treatment. 
#' 
#' logit.active= the logit of the probability of an outcome in the active treatment.
#' 
#' benefit= treatment benefit in log odds ratio. 
#' 
#' py=the probability of the outcome for each patient, under the treatment actually administered.
#' 
#' logit.py= the logit of py.
#' 
#' y.observed= the observed outcome
#' @examples 
#' dat1=simbinary(100)$dat
#' head(dat1)
#' @export
simbinary=function(Npat=100){
  ### simulate covariates
  library(MASS)
  expit=function(x){ exp(x)/(1 + exp(x))}
logit=function(x){l=log(x/(1-x)); return(l)}
  # observed true predictors
  x1x2<-mvrnorm(n = Npat, c(0,0), Sigma=matrix(c(1,0.2,0.2,1),2,2))
  simdat <- data.frame(x1=x1x2[,1])
  simdat$x2<-x1x2[,2] 
  x3.lin <-rbinom(Npat, 1, prob = 0.2)
  simdat$x3<-rbinom(Npat, 1, prob = expit(x3.lin) )
  
  #unobserved true predictors
  simdat$x4<-rnorm(Npat, 0, 1 )
  
  ## observed nuisance parameters
  simdat$x5<-rnorm(Npat, 0, 1 )
  simdat$x6<-rnorm(Npat, 0, 1 )

  
  ### simulate treatment assignment
  pt <- 0.5 ## randomized
  simdat$t <- rbinom(Npat, 1, prob = pt)
  ### outcome at control treatment
  simdat$logit.control<- with(simdat, -2+0.2*x1+1*x2+0.6*x3+0.3*x4+
                            0.1*x1^2+0.2*x1*x4+0.2*x4*x4+
                            rnorm(Npat, 0, 0.1))
  ### treatment benefit
  simdat$benefit<-with(simdat, -0.5-0.2*x1-0.1*x2+0.2*x3+0.2*x4+0.5*x1*x4+
                         0.1*x4*x2+0.1*x2*x3+0.1*x3*x4+
                         rnorm(Npat, 0, 0.05))
  simdat$logit.active=with(simdat,logit.control+benefit)
  simdat$logit.py=with(simdat,logit.control+t*benefit)
  simdat$py<-expit(simdat$logit.py)
  simdat$y.observed=rbinom(Npat, 1, prob=simdat$py)
  return(list(dat=simdat))
}


#' Expit
#'
#' Calculates the expit of a real number
#' @param x A real number
#' @return exp(x)/(1+exp(x))
#' @examples 
#' expit(2.3)
#' @export
expit=function(x){ exp(x)/(1 + exp(x))}


#' Logit
#'
#' Calculates the expit of a real number between 0 and 1
#' @param x A real number between 0 and 1
#' @return log(x/(1-x))
#' @examples 
#' logit(0.2)
#' @export
logit=function(x){l=log(x/(1-x)); return(l)}