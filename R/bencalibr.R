#' Plotting calibration for benefit of a prediction model
#'
#' This function produces plot to illustrate the calibration for benefit 
#' for a prediction model. The data are split in a number of groups according
#' to the predicted benefit, and within each group the functions estimates
#' the observed treatment benefit and compares it with the predicted one 
#' according to the prediction model
#' @param Ngroups The number of groups to split the data.
#' @param y.observed The observed outcome.
#' @param treat A vector with the treatment assignment. This must be 0 (for control treatment)
#' or 1 (for active treatment).
#' @param predicted.treat.0 A vector with the model predictions for each patient, under the control treatment.
#' For the case of a binary outcome this should be in the logit scale.
#' @param predicted.treat.1 A vector with the model predictions for each patient, under the active treatment.
#' For the case of a binary outcome this should be in the logit scale.
#' @param type The type of the outcome, "binary" or "continuous".
#' @param smoothing.function The method used to smooth the calibration line. Can 
#' be "lm", "glm", "gam", "loess", "rlm". More details can be found in https://ggplot2.tidyverse.org/reference/geom_smooth.html.
#' @param measure For binary outcomes only. Can be risk difference ("RD") or 
#' log odds ratios ("logor")
#' @param axis.limits Sets the limits of the graph. This can be omitted.
#' @return The calibration plot
#' @examples 
#' # continuous outcome 
#' dat1=simcont(200)$dat
#' head(dat1)
#' lm1=lm(y.observed~(x1+x2+x3)*t, data=dat1)
#' dat.t0=dat1; dat.t0$t=0 
#' dat.t1=dat1; dat.t1$t=1
#' dat1$predict.treat.1=predict(lm1, newdata = dat.t1) # predictions in treatment
#' dat1$predict.treat.0=predict(lm1, newdata = dat.t0) # predicions in control

#' bencalibr(data=dat1, Ngroups=10, y.observed, predicted.treat.1=predict.treat.1,
#'           predicted.treat.0=predict.treat.0, type="continuous", treat=t, 
#'           smoothing.function = "lm", axis.limits = c(-1, 1.3))
#' # binary outcome 
#' dat2=simbinary(300)$dat
#' head(dat2)
#' glm1=glm(y.observed~(x1+x2+x3)*t, data=dat2, family = binomial(link = "logit"))
#' dat2.t0=dat2; dat2.t0$t=0 
#' dat2.t1=dat2; dat2.t1$t=1
#' dat2$predict.treat.1=predict(glm1, newdata = dat2.t1) # predictions in treatment
#' dat2$predict.treat.0=predict(glm1, newdata = dat2.t0) # predicions in control
#' bencalibr(data=dat2, Ngroups=6, y.observed, predicted.treat.1=predict.treat.1,
#'           predicted.treat.0=predict.treat.0, type="binary", treat=t, 
#'           smoothing.function = "lm", measure="logor")
#' bencalibr(data=dat2, Ngroups=6, y.observed, predicted.treat.1=predict.treat.1,
#'           predicted.treat.0=predict.treat.0, type="binary", treat=t, 
#'           smoothing.function = "lm", measure="RD")
#' @export
bencalibr=function(data=NULL, Ngroups=5, y.observed, treat, 
                   predicted.treat.0, predicted.treat.1 ,type="continuous", 
                   smoothing.function="rlm", measure = "RD", axis.limits=NULL){
  ##smoothing function can be: lm, glm, gam, loess, rlm
  ##measure can be RD or logor
  library(ggplot2)
  se <- function(x) sd(x)/sqrt(length(x))
  nulldata=is.null(data)
  
  ### continuous ----
  if(type=="continuous"){
    
    if(nulldata){
      data<-data.frame(y.observed, treat,predicted.treat.0, predicted.treat.1 ) }
    
    if(!nulldata){
      arguments <- as.list(match.call())
      y.observed = eval(arguments$y.observed, data)
      treat=eval(arguments$treat, data)
      predicted.treat.0=eval(arguments$predicted.treat.0, data)
      predicted.treat.1=eval(arguments$predicted.treat.1, data)
      data<-data.frame(y.observed, treat=treat,predicted.treat.0, predicted.treat.1 ) 
      data$predicted.benefit<-data$predicted.treat.1-data$predicted.treat.0}
    
    data$predicted.benefit<-data$predicted.treat.1-data$predicted.treat.0
    data$predicted.y<- data$predicted.treat.0*(data$treat==0)+ data$predicted.treat.1*(data$treat==1)
    d1<-quantile(data$predicted.benefit, probs = seq(0, 1, 1/Ngroups))
    
    g1<-list()
    for (i in 1:Ngroups){g1[[i]]<-data[data$predicted.benefit>=d1[i]&data$predicted.benefit<d1[i+1],]}
    
    predicted=c()
    observed=c()
    errors.predicted=c()
    errors.observed=c()
    for (i in 1:Ngroups){
      predicted=c(predicted,mean(g1[[i]]$predicted.benefit))
      observed=c(observed, 
                 mean(g1[[i]]$y.observed[g1[[i]]$treat==1])-mean(g1[[i]]$y.observed[g1[[i]]$treat==0])
      )
      
      errors.predicted=c(errors.predicted,sd(g1[[i]]$predicted.benefit))
      errors.observed=c(errors.observed, 
                        sqrt(se(g1[[i]]$y.observed[g1[[i]]$t==1])^2+
                               se(g1[[i]]$y.observed[g1[[i]]$t==0])^2)
      )  }
    
    dat1=data.frame("pred"=predicted, "obs"=observed, "se.pred"=errors.predicted, "se.obs"=errors.observed)
    
    
    p1= suppressMessages(ggplot(dat1, aes(x=pred, y=obs)) +
                           geom_point(size=3, shape=20)+
                           labs(x="Predicted benefit", y="Observed benefit")+ 
                           
                           geom_abline(intercept = 0, slope = 1, color="black", 
                                       linetype="dashed", size=0.5)+
                           geom_errorbar(aes(ymin=obs-se.obs, ymax=obs+se.obs), width=.1) +
                           geom_errorbarh(aes(xmin=pred-se.pred, xmax=pred+se.pred))+
                           geom_smooth(method=smoothing.function, colour="blue",size=0.5)+ theme(aspect.ratio=1)
    )
    r1<-max((suppressMessages(layer_scales(p1))$x$range$range))
    r2<-min((suppressMessages(layer_scales(p1)$x$range$range)))
    s1<-max((suppressMessages(layer_scales(p1)$y$range$range)))
    s2<-min((suppressMessages(layer_scales(p1)$y$range$range)))
    t1<-(max(r1,s1))
    t2<-(min(r2,s2))
    if(is.null(axis.limits)){t11=t1; t12=t2}
    if(!is.null(axis.limits)){t11=axis.limits[2]; t12=axis.limits[1]}
    p1<- p1+coord_equal(xlim=c(t12,t11),ylim=c(t12,t11))
    
    return(suppressMessages(print(p1)))
    cat(" Type of outcome: ", type, "\n",
        "Ngroups: ", Ngroups,"\n",
        "Smoothing function: ",smoothing.function,"\n"    )
    
  }
  #### binary ----
  if (type=="binary"){
    logit<-function(x){l<-log(x/(1-x)); return(l)}
    if(nulldata){data<-data.frame(y.observed, treat,predicted.treat.0, predicted.treat.1 ) }
    if(!nulldata){
      arguments <- as.list(match.call())
      y.observed = eval(arguments$y.observed, data)
      treat=eval(arguments$treat, data)
      predicted.treat.0=eval(arguments$predicted.treat.0, data)
      predicted.treat.1=eval(arguments$predicted.treat.1, data)
      data<-data.frame(y.observed, treat, predicted.treat.0, predicted.treat.1) 
    }
    #### logor ----
    if (measure=="logor"){
      
      data$benefit<-data$predicted.treat.1-data$predicted.treat.0
      d1<-quantile(data$benefit, probs <- seq(0, 1, 1/Ngroups))
      g1<-list()
      for (i in 1:Ngroups){g1[[i]]<-data[data$benefit>=d1[i]&data$benefit<d1[i+1],]}
      predicted<-c()
      observed<-c()
      errors.predicted<-c()
      errors.observed<-c()
      for (i in 1:Ngroups){
        predicted=c(predicted,mean(g1[[i]]$benefit))
        observed=c(observed, 
                   logit(mean(g1[[i]]$y.observed[g1[[i]]$treat==1]))-
                     logit(mean(g1[[i]]$y.observed[g1[[i]]$treat==0])))
        errors.predicted<-c(errors.predicted,sd(g1[[i]]$benefit))
        errors.observed<-c(errors.observed, 
                           sqrt(
                             1/sum(g1[[i]]$y.observed[g1[[i]]$treat==1])
                             +1/sum(g1[[i]]$y.observed[g1[[i]]$treat==1]==0)
                             +1/sum(g1[[i]]$y.observed[g1[[i]]$treat==0])
                             +1/sum(g1[[i]]$y.observed[g1[[i]]$treat==0]==0) ))
      }
      
      dat1<-data.frame("pred"<-predicted, "obs"<-observed, 
                       "se.pred"<-errors.predicted, "se.obs"<-errors.observed)
      
      p1<- ggplot(dat1, aes(x<-pred, y<-obs)) +
        geom_point(size=3, shape=20)+
        labs(x="Predicted benefit", y="Observed benefit")+ 
        
        geom_abline(intercept = 0, slope = 1, color="black", 
                    linetype="dashed", size=0.5)+
        geom_errorbar(aes(ymin=obs-se.obs, ymax=obs+se.obs), width=.1) +
        geom_errorbarh(aes(xmin=pred-se.pred, xmax=pred+se.pred)) +#+ coord_cartesian( xlim = c(-1, 1),  ylim = c(-1, 1))
        geom_smooth(method="lm", colour="blue",size=0.5)+ theme(aspect.ratio=1)+
        ggtitle("Calibration for benefit (measure: log odds ratios)")
      
      r1<-max((suppressMessages(layer_scales(p1))$x$range$range))
      r2<-min((suppressMessages(layer_scales(p1)$x$range$range)))
      s1<-max((suppressMessages(layer_scales(p1)$y$range$range)))
      s2<-min((suppressMessages(layer_scales(p1)$y$range$range)))
      
      t1<-(max(r1,s1))
      t2<-(min(r2,s2))
      if(is.null(axis.limits)){t11=t1; t12=t2}
      if(!is.null(axis.limits)){t11=axis.limits[2]; t12=axis.limits[1]}
      p1<- p1+coord_equal(xlim=c(t12,t11),ylim=c(t12,t11))
      
      cat(" Type of outcome: ", type,"\n", 
          "Ngroups: ", Ngroups,"\n",
          "Smoothing.function: ",smoothing.function,"\n", 
          "Type of measure:", measure, "\n"
      )
      
      return(suppressMessages(print(p1)))
      
    }
    
    #### rd ----
    if (measure=="RD"){
      expit<-function(x){ exp(x)/(1 + exp(x))}
      data$benefit<-expit(data$predicted.treat.1)-expit(data$predicted.treat.0)
      d1<-quantile(data$benefit, probs <- seq(0, 1, 1/Ngroups))
      g1<-list()
      for (i in 1:Ngroups){g1[[i]]<-data[data$benefit>=d1[i]&data$benefit<d1[i+1],]}
      predicted<-c()
      observed<-c()
      errors.predicted<-c()
      errors.observed<-c()
      for (i in 1:Ngroups){
        predicted=c(predicted,mean(g1[[i]]$benefit))
        
        observed=c(observed, 
                   (mean(g1[[i]]$y.observed[g1[[i]]$t==1]))-
                     (mean(g1[[i]]$y.observed[g1[[i]]$t==0])))
        
        errors.predicted<-c(errors.predicted,sd(expit(g1[[i]]$benefit)))
        errors.observed<-c(errors.observed, 
                           sqrt(
                             sum(g1[[i]]$y.observed[g1[[i]]$t==1])*
                               sum(g1[[i]]$y.observed[g1[[i]]$t==1]==0)/
                               length(g1[[i]]$y.observed[g1[[i]]$t==1])^3+
                               
                               sum(g1[[i]]$y.observed[g1[[i]]$t==0])*
                               sum(g1[[i]]$y.observed[g1[[i]]$t==0]==0)/
                               length(g1[[i]]$y.observed[g1[[i]]$t==0])^3
                           ))
      }
      
      dat1<-data.frame("pred"<-predicted, "obs"<-observed, "se.pred"<-errors.predicted, "se.obs"<-errors.observed)
      
      p2<- ggplot(dat1, aes(x<-pred, y<-obs)) +
        geom_point(size=3, shape=20)+
        labs(x="Predicted benefit", y="Observed benefit")+ 
        geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
        geom_errorbar(aes(ymin=obs-se.obs, ymax=obs+se.obs), width=.005) +
        geom_errorbarh(aes(xmin=pred-se.pred, xmax=pred+se.pred), height=.005) +
        geom_smooth(method="lm", colour="blue",size=0.5)+ theme(aspect.ratio=1)+
        ggtitle("Calibration for benefit (measure: risk difference)")
      
      r1<-max((suppressMessages(layer_scales(p2))$x$range$range))
      r2<-min((suppressMessages(layer_scales(p2)$x$range$range)))
      s1<-max((suppressMessages(layer_scales(p2)$y$range$range)))
      s2<-min((suppressMessages(layer_scales(p2)$y$range$range)))
      
      t1<-(max(r1,s1))
      t2<-(min(r2,s2))
      if(is.null(axis.limits)){t11=t1; t12=t2}
      if(!is.null(axis.limits)){t11=axis.limits[2]; t12=axis.limits[1]}
      p2<- p2+coord_equal(xlim=c(t12,t11),ylim=c(t12,t11))
      cat(" Type of outcome: ", type,"\n", 
          "Ngroups: ", Ngroups,"\n",
          "Smoothing.function: ",smoothing.function,"\n",
          "Type of measure:", measure, "\n"
      )
      
      return(suppressMessages(print(p2)))
    }
    
  }
}


