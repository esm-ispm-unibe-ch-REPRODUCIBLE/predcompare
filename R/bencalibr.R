#' Calibration for benefit plots
#'
#' Produces plot
#' @param data A dataframe
#' @return Aplot
#' @examples 
#' temp1 <- F_to_C(50);
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export
bencalibr=function(data=NULL, Ngroups=5, y.observed, treat, 
                   predicted.treat.0, predicted.treat.1 ,type="continuous", 
                   smoothing.function="rlm", measure = "RD"){
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
    p1<-suppressMessages(p1+coord_equal(xlim=c(t2,t1),ylim=c(t2,t1)))
  
    
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
    p1<- p1+coord_equal(xlim=c(t2,t1),ylim=c(t2,t1))
   
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
      p2<- p2+coord_equal(xlim=c(t2,t1),ylim=c(t2,t1))
      cat(" Type of outcome: ", type,"\n", 
          "Ngroups: ", Ngroups,"\n",
          "Smoothing.function: ",smoothing.function,"\n",
          "Type of measure:", measure, "\n"
      )
      
      return(suppressMessages(print(p2)))
    }
    
    }
}

