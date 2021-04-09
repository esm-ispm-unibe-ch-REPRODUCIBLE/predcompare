#' Calculating measures for calibration for benefit for a prediction model
#'
#' This function calculates a series of measures to assess the calibration for benefit of a prediction model.
#' For a continuous outcome it uses a one-to-one matching of patients based on covariates or 
#' predicted treatment benefit. It also uses a matching in groups using the k-means algorithm. It
#' also fits a regression line for benefit. For a binary outcome it also calculates the c-for-benefit. 
#' @param repeats The number of repetitions for the algorithm.
#' @param Ngroups The number of groups to split the data.
#' @param X A dataframe with patient covariates.
#' @param Y The observed outcome.
#' @param treat A vector with the treatment assignment. This must be 0 (for control treatment)
#' or 1 (for active treatment).
#' @param predicted.treat.0 A vector with the model predictions for each patient, under the control treatment.
#' For the case of a binary outcome this should be in the logit scale.
#' @param predicted.treat.1 A vector with the model predictions for each patient, under the active treatment.
#' For the case of a binary outcome this should be in the logit scale.
#' @param type The type of the outcome, "binary" or "continuous".
#' @param measure For binary outcomes only. Can be risk difference ("RD") or 
#' log odds ratios ("logor")
#' @return A table
#' @examples 
#' # continuous outcome 
#' dat1=simcont(500)$dat
#' lm1=lm(y.observed~(x1+x2+x3)*t, data=dat1)
#' dat.t0=dat1; dat.t0$t=0 # a dataset 
#' dat.t1=dat1; dat.t1$t=1
#' dat1$predict.treat.1=predict(lm1, newdata = dat.t1) # predictions in treatment
#' dat1$predict.treat.0=predict(lm1, newdata = dat.t0) # predicions in control
#' 
#' predcompare(repeats=20, Ngroups=c(5:10), X=dat1[,c("x1", "x2","x3")], 
#'             Y=dat1$y.observed, 
#'             predicted.treat.1 = dat1$predict.treat.1,
#'             predicted.treat.0 = dat1$predict.treat.0,
#'             treat=dat1$t, type="continuous")
#'             
#' # binary outcome 
#' dat2=simbinary(800)$dat
#' head(dat2)
#' glm1=glm(y.observed~(x1+x2+x3)*t, data=dat2, family = binomial(link = "logit"))
#' dat2.t0=dat2; dat2.t0$t=0 
#' dat2.t1=dat2; dat2.t1$t=1
#' dat2$predict.treat.1=predict(glm1, newdata = dat2.t1) # predictions in treatment
#' dat2$predict.treat.0=predict(glm1, newdata = dat2.t0) # predicions in control
#' 
#' predcompare(repeats=20, Ngroups=c(5:10), X=dat2[,c("x1", "x2","x3")], 
#'             Y=dat2$y.observed, 
#'             predicted.treat.1 = dat2$predict.treat.1,
#'             predicted.treat.0 = dat2$predict.treat.0,
#'             treat=dat2$t, type="binary", measure="RD")
#'             
#' predcompare(repeats=20, Ngroups=c(5:10), X=dat2[,c("x1", "x2","x3")], 
#'             Y=dat2$y.observed, 
#'             predicted.treat.1 = dat2$predict.treat.1,
#'             predicted.treat.0 = dat2$predict.treat.0,
#'             treat=dat2$t, type="binary", measure="logor")
#' @export
predcompare=function(repeats=50, ### number of repeated uses of the algorithm
                     Ngroups=10, ### number of groups for the k-means analysis
                     X, ### a data frame of the covariates
                     treat, ### a vector of treatment assignment. Should be 1 or 0  
                     Y,  ### a vector of the outcomes. For binary this should be 1 or 0
                     predicted.treat.1, ### a vector of the predicted outcome in treatment 1. 
                     # for binary outcomes this should be a log odds
                     predicted.treat.0, ### a vector of the predicted outcome in treatment 0. 
                     # for binary outcomes this should be a log odds
                     type="continuous", ### either continuous or binary
                     measure="RD" ### for binary outcomes only; either "RD" (risk difference) 
                     ### or "logit" is allowed
){
  library(Matching)
  library(Hmisc)
  Ngroups.length=length(Ngroups)
  
  ################ continuous data --------------
  if(type=="continuous"){
    dat1<-cbind(X,"Y"=c(Y),"t"=treat,"benefit"=predicted.treat.1-predicted.treat.0,
                "predicted.t1"=predicted.treat.1, "predicted.t0"=predicted.treat.0)
    # k-means
    bias.m1=mse.m1=coeff.m11=coeff.m12=rsq.m1=c()
    results.kmeans<-list()
    for(n in 1:Ngroups.length){
      
      
      for (k in 1:repeats)
      {
        groups.kmean<-kmeans(x=dat1[,colnames(X)], centers = Ngroups[n], iter.max = 50)
        dat1$group<-groups.kmean$cluster
        ngroups<-as.vector(table(dat1$group))
        observed.groups1<-c()
        observed.groups0<-c()
        estimated.groups.m1<-c()
        dall<-NULL
        for( i in 1:Ngroups[n]){
          # observed
          observed.groups1=c(observed.groups1,(mean(dat1$Y[dat1$group==i & dat1$t==1])))
          observed.groups0=c(observed.groups0,(mean(dat1$Y[dat1$group==i & dat1$t==0])))
          #predicted
          
          estimated.groups.m1=c(estimated.groups.m1,mean(dat1$predicted.t1[dat1$group==i&dat1$t==1])-mean(dat1$predicted.t0[dat1$group==i&dat1$t==0])) 
          
        } 
        # collate resuts
        dall<-data.frame(observed.groups0, observed.groups1, estimated.groups.m1,ngroups)
        dall<-dall[complete.cases(dall),]
        dall$observed.groups=dall$observed.groups1-dall$observed.groups0
        reg.m1<-summary(lm(dall$observed.groups~dall$estimated.groups.m1, weights= dall$ngroups))
        coeff.m11<-c( coeff.m11, coef(reg.m1)[1])
        coeff.m12<-c( coeff.m12,coef(reg.m1)[2])
        rsq.m1<-c(rsq.m1, reg.m1$r.squared)
        bias.m1<-c(bias.m1, mean(dall$observed.groups-dall$estimated.groups.m1))
        mse.m1<-c(mse.m1, mean((dall$observed.groups-dall$estimated.groups.m1)^2)) }
      # summary
      
      results.kmeans[[n]]<-data.frame(bias= round(median(bias.m1),digits=2),
                                      mse= round(median(mse.m1),digits=2),
                                      a0= round(median(coeff.m11),digits=2),
                                      a1= round(median(coeff.m12),digits=2),
                                      R2= round(median(rsq.m1),digits=2)  )        
    }
    
    # match by X ---
    bias.m1=mse.m1=m1.a0=m1.a1=m1.R2=c()
    for(k in 1:repeats){
      g1=dat1[dat1$t==1,]
      g0=dat1[dat1$t==0,]
      n1=min(length(g1$t) ,length(g0$t))
      g1=g1[sample(length(g1$t),n1),]
      g0=g0[sample(length(g0$t),n1),]
      dataall=rbind(g1,g0)
      covariates<-dataall[,colnames(X)]
      
      matched=Match(Tr=dataall$t, X=covariates)
      
      data.compare=data.frame(tr.obs=dataall$Y[matched$index.treated], ctr.obs=dataall$Y[matched$index.control], 
                              benefit.m1=(dataall$benefit[matched$index.treated]+dataall$benefit[matched$index.control])/2)
      
      data.compare$obs.benefit=with(data.compare, tr.obs-ctr.obs)
      bias.m1=c(bias.m1, mean(data.compare$obs.benefit-data.compare$benefit.m1))
      data.compare$obs.benefit=with(data.compare, tr.obs-ctr.obs)
      mse.m1=c(mse.m1, mean((data.compare$obs.benefit-data.compare$benefit.m1)^2))
      reg.m1=lm(data.compare$obs.benefit~data.compare$benefit.m1)
      m1.a0=c(m1.a0, coef(reg.m1)[1]); m1.a1=c(m1.a1, coef(reg.m1)[2]);  m1.R2=c(m1.R2, summary(reg.m1)$r.squared)
    }
    
    # summary
    results.after.matching=data.frame(bias=round(median(bias.m1),digits=2),
                                      mse=round(median(mse.m1),digits=2),
                                      a0= round(median(m1.a0),digits=2),
                                      a1= round(median(m1.a1),digits=2),
                                      R2= round(median(m1.R2),digits=2)                      )
    
    ### match by benefit----
    dataall=g0=g1=n1=covariates=data.compare=matched=reg.m1=NULL
    bias.m1=mse.m1=m1.a0=m1.a1=m1.R2=c()
    for(k in 1:repeats){
      g1=dat1[dat1$t==1,]
      g0=dat1[dat1$t==0,]
      n1=min(length(g1$t) ,length(g0$t))
      g1=g1[sample(length(g1$t),n1),]
      g0=g0[sample(length(g0$t),n1),]
      dataall=rbind(g1,g0)
      matched1=Match(Tr=dataall$t, X=dataall$benefit)
      
      data.compare=data.frame(tr.obs=dataall$Y[matched1$index.treated], ctr.obs=dataall$Y[matched1$index.control], 
                              benefit.m1=(dataall$benefit[matched1$index.treated]+dataall$benefit[matched1$index.control])/2)
      
      data.compare$obs.benefit=with(data.compare, tr.obs-ctr.obs)
      bias.m1=c(bias.m1, mean(data.compare$obs.benefit-data.compare$benefit.m1))
      data.compare$obs.benefit=with(data.compare, tr.obs-ctr.obs)
      mse.m1=c(mse.m1, mean((data.compare$obs.benefit-data.compare$benefit.m1)^2))
      reg.m1=lm(data.compare$obs.benefit~data.compare$benefit.m1)
      m1.a0=c(m1.a0, coef(reg.m1)[1]); m1.a1=c(m1.a1, coef(reg.m1)[2]);  m1.R2=c(m1.R2, summary(reg.m1)$r.squared)
    }
    ## summary
    results.after.matching.benefit=data.frame(bias=round(median(bias.m1),digits=2),
                                              mse=round(median(mse.m1),digits=2),
                                              a0= round(median(m1.a0),digits=2),
                                              a1= round(median(m1.a1),digits=2),
                                              R2= round(median(m1.R2),digits=2)                      )
    
    results=data.frame(Method=c("match by covariates", "match by benefit", paste("kmeans N=",Ngroups[1], sep="")), 
                       Bias=c( results.after.matching$bias ,results.after.matching.benefit$bias, results.kmeans[[1]]$bias)
                       ,MSE=c( results.after.matching$mse ,results.after.matching.benefit$mse, results.kmeans[[1]]$mse)
                       ,a0=c( results.after.matching$a0 ,results.after.matching.benefit$a0, results.kmeans[[1]]$a0)
                       ,a1=c( results.after.matching$a1 ,results.after.matching.benefit$a1, results.kmeans[[1]]$a1)
                       ,R2=c( results.after.matching$R2 ,results.after.matching.benefit$R2, results.kmeans[[1]]$R2) )
    
    if(Ngroups.length>1){
      
      
      for (n in 2:Ngroups.length){
        d1<-data.frame(Method=c(paste("kmeans N=",Ngroups[n], sep="")), 
                       Bias=c(results.kmeans[[n]]$bias)
                       ,MSE=c( results.kmeans[[n]]$mse)
                       ,a0=c(  results.kmeans[[n]]$a0)
                       ,a1=c(results.kmeans[[n]]$a1)
                       ,R2=c( results.kmeans[[n]]$R2))
        results=rbind(results, d1)}
      
      
    }  
    dat1$predicted.treat.0<-predicted.treat.0
    lm1<-lm(dat1$Y~-1+dat1$predicted.treat.0+dat1$t:dat1$benefit)
    sum.lm1<-summary(lm1)
    c1=format(round(coef(lm1)[1], 2), nsmall = 2)
    c2=format(round(coef(lm1)[2], 2), nsmall = 2)
    R2=format(round(sum.lm1$r.squared, 2), nsmall = 2)
    
    results.reg=data.frame("Intercept"=c1, "Slope"=c2, "R2"=R2)
    rownames(results.reg)=c("Estimates") 
    results.all=list("matching.methods"=results, "regression.for.benefit"=results.reg)
    
    cat(" Type of outcome: ", type, "\n",
        "Repeats: ",repeats,"\n"    )
    
  }
  
  
  ################ binary data  --------------
  
  if(type=="binary"){
    library(rsq)
    cat(" Type of outcome: binary","\n", 
        "Repeats: ", repeats, "\n", 
        "Measure: ", measure, "\n")   
    ### c for benefit --------
    dat1<-cbind(X,"Y"=c(Y),"t"=treat, "benefit"=c(predicted.treat.1-predicted.treat.0))
    c.by.benefit=c(); se.c.by.benefit=c()
    c.by.covariates=c(); se.c.by.covariates=c()
    n.0 <- sum(dat1$t==0)
    n.1 <- sum(dat1$t==1)
    n.10=min(n.0, n.1)
    for (i in 1:repeats){
      rows.in<-c(which(dat1$t==1)[sample(n.1, n.10)], which(dat1$t==0)[sample(n.0, n.10)])
      dat2=dat1[rows.in[order(rows.in)],]
      X2=X[rows.in[order(rows.in)],]
      ## match by benefit
      ind.0 <- which(dat2$t==0)
      order.0 <- order(dat2$benefit[ind.0])
      ind.0 <- ind.0[order.0]
      
      ind.1 <- which(dat2$t==1)
      order.1 <- order(dat2$benefit[ind.1])
      ind.1 <- ind.1[order.1]
      
      pred.ben.0 <- dat2$benefit[ind.0]
      pred.ben.1 <- dat2$benefit[ind.1]
      pred.ben.avg <- (pred.ben.1+pred.ben.0)/2
      obs.out.0 <- dat2$Y[ind.0]
      obs.out.1<- dat2$Y[ind.1]
      obs.ben <- obs.out.1-obs.out.0
      
      cindex <- rcorr.cens(pred.ben.avg, obs.ben)
      c.by.benefit <- c(c.by.benefit, cindex["C Index"][[1]])
      se.c.by.benefit<-c(se.c.by.benefit, cindex["S.D."][[1]]/2	)
      
      ## match by covariates
      rr1 <- Match(Tr=dat2$t, X=X2, M=1,ties=F,replace=FALSE)
      ind.0 <- rr1$index.control
      ind.1 <- rr1$index.treated
      ### Calculation of predicted and observed benefit in matched pairs
      pred.ben.0 <- dat2$benefit[ind.0]
      pred.ben.1 <- dat2$benefit[ind.1]
      pred.ben.avg <- (pred.ben.1+pred.ben.0)/2
      obs.out.0 <- dat2$Y[ind.0]
      obs.out.1 <- dat2$Y[ind.1]
      obs.ben <- obs.out.1-obs.out.0
      # Benefit c-statistic
      cindex2 <- rcorr.cens(pred.ben.avg, obs.ben)
      c.by.covariates <- c(c.by.covariates, cindex2["C Index"][[1]])
      se.c.by.covariates<-c(se.c.by.covariates, cindex2["S.D."][[1]]/2	)
      
    }
    c.by.covariates1=format(round(mean(c.by.covariates), 2), nsmall = 2)
    se.c.by.covariates1=format(round(mean(se.c.by.covariates), 2), nsmall = 2)
    
    c.by.benefit1=format(round(mean(c.by.benefit), 2), nsmall = 2)
    se.c.by.benefit1=format(round(mean(se.c.by.benefit), 2), nsmall = 2)
    results2=data.frame("c.for.benefit"=c(c.by.benefit1, c.by.covariates1), 
                        "SE"=c(se.c.by.benefit1,se.c.by.covariates1) )
    rownames(results2)=c("Matched by benefit", "Matched by covariates")
    
    ###
    
    logit<-function(x){l=log(x/(1-x)); return(l)}
    ### logor ---- 
    bias.m=mse.m=coeff.1=coeff.2=rsq.m=c()
    results.kmeans.bin<-list()
    if (measure=="logor"){
      dat1$predicted.t1=expit(predicted.treat.1)
      dat1$predicted.t0=expit(predicted.treat.0)
      for(n in 1:Ngroups.length){
        for (k in 1:repeats)
        {
          groups.kmean<-kmeans(x=dat1[,colnames(X)], centers = Ngroups[n], iter.max = 50)
          dat1$group<-groups.kmean$cluster
          ngroups<-as.vector(table(dat1$group))
          observed.groups1<-c()
          observed.groups0<-c()
          estimated.groups<-c()
          dall=NULL
          for( i in 1:Ngroups[n]){
            # observed
            observed.groups1=c(observed.groups1, 
                               (logit(mean(dat1$Y[dat1$group==i & dat1$t==1]))))
            observed.groups0=c(observed.groups0, 
                               (logit(mean(dat1$Y[dat1$group==i & dat1$t==0]))))
            # estimated
            estimated.groups=c(estimated.groups,
                               logit(mean(dat1$predicted.t1[dat1$group==i & dat1$t==1]))-
                                 logit(mean(dat1$predicted.t0[dat1$group==i & dat1$t==0]))) 
          } 
          # collate resuts
          dall=data.frame(observed.groups0, observed.groups1, estimated.groups) 
          is.na(dall$observed.groups0[is.infinite(dall$observed.groups0)])=TRUE
          is.na(dall$observed.groups1[is.infinite(dall$observed.groups1)])=TRUE
          dall=dall[complete.cases(dall),]
          dall$observed.groups=dall$observed.groups1-dall$observed.groups0
          reg=summary(lm(dall$observed.groups~dall$estimated.groups))
          coeff.1=c( coeff.1, coef(reg)[1])
          coeff.2=c( coeff.2,coef(reg)[2])
          rsq.m=c(rsq.m, reg$r.squared)
          bias.m=c(bias.m, mean(dall$observed.groups-dall$estimated.groups))
          mse.m=c(mse.m, mean((dall$observed.groups-dall$estimated.groups)^2))
        }
        # summary
        results.kmeans.bin[[n]]<-data.frame(bias= round(median(bias.m),digits=2),
                                            mse= round(median(mse.m),digits=2),
                                            a0= round(median(coeff.1),digits=2),
                                            a1= round(median(coeff.2),digits=2),
                                            R2= round(median(rsq.m),digits=2)  )  
      }
      
      results=data.frame(Method=c(paste("kmeans N=",Ngroups[1], sep="")), 
                         Bias=c( results.kmeans.bin[[1]]$bias)
                         ,MSE=c(results.kmeans.bin[[1]]$mse)
                         ,a0=c(  results.kmeans.bin[[1]]$a0)
                         ,a1=c( results.kmeans.bin[[1]]$a1)
                         ,R2=c( results.kmeans.bin[[1]]$R2))
      
      if(Ngroups.length>1){
        for (n in 2:Ngroups.length){
          d1<-data.frame(Method=c(paste("kmeans N=",Ngroups[n], sep="")), 
                         Bias=c(results.kmeans.bin[[n]]$bias)
                         ,MSE=c( results.kmeans.bin[[n]]$mse)
                         ,a0=c(  results.kmeans.bin[[n]]$a0)
                         ,a1=c(results.kmeans.bin[[n]]$a1)
                         ,R2=c( results.kmeans.bin[[n]]$R2))
          results=rbind(results, d1)
          
        }
      } 
    }
    
    
    ### RD ---- 
    if (measure=="RD"){
      expit<-function(x){ exp(x)/(1 + exp(x))}
      dat1<-cbind(X,"Y"=Y,"t"=treat, "predicted.t1"=expit(predicted.treat.1),
                  "predicted.t0"=expit(predicted.treat.0))
      dat1$benefit<-(dat1$predicted.t1)-(dat1$predicted.t0)
      
      for(n in 1:Ngroups.length){
        for (k in 1:repeats)
        {
          groups.kmean<-kmeans(x=dat1[,colnames(X)], centers = Ngroups[n], iter.max = 50)
          dat1$group<-groups.kmean$cluster
          ngroups<-as.vector(table(dat1$group))
          observed.groups1<-c()
          observed.groups0<-c()
          estimated.groups<-c()
          dall=NULL
          for( i in 1:Ngroups[n]){
            # observed
            observed.groups1=c(observed.groups1, 
                               (mean(dat1$Y[dat1$group==i & dat1$t==1])))
            
            observed.groups0=c(observed.groups0, 
                               (mean(dat1$Y[dat1$group==i & dat1$t==0])))
            # estimated
            estimated.groups=c(estimated.groups,mean(dat1$predicted.t1[dat1$group==i&dat1$t==1])-mean(dat1$predicted.t0[dat1$group==i&dat1$t==0])) 
            
          } 
          # collate resuts
          dall=data.frame(observed.groups0, observed.groups1, estimated.groups) 
          is.na(dall$observed.groups0[is.infinite(dall$observed.groups0)])=TRUE
          is.na(dall$observed.groups1[is.infinite(dall$observed.groups1)])=TRUE
          dall=dall[complete.cases(dall),]
          dall$observed.groups=dall$observed.groups1-dall$observed.groups0
          reg=summary(lm(dall$observed.groups~dall$estimated.groups))
          coeff.1=c( coeff.1, coef(reg)[1])
          coeff.2=c( coeff.2,coef(reg)[2])
          rsq.m=c(rsq.m, reg$r.squared)
          bias.m=c(bias.m, mean(dall$observed.groups-dall$estimated.groups))
          mse.m=c(mse.m, mean((dall$observed.groups-dall$estimated.groups)^2))
        }
        # summary
        results.kmeans.bin[[n]]<-data.frame(bias= round(median(bias.m),digits=3),
                                            mse= round(median(mse.m),digits=3),
                                            a0= round(median(coeff.1),digits=2),
                                            a1= round(median(coeff.2),digits=2),
                                            R2= round(median(rsq.m),digits=2)  )  
      }
      
      results=data.frame(Method=c(paste("kmeans N=",Ngroups[1], sep="")), 
                         Bias=c( results.kmeans.bin[[1]]$bias)
                         ,MSE=c(results.kmeans.bin[[1]]$mse)
                         ,a0=c(  results.kmeans.bin[[1]]$a0)
                         ,a1=c( results.kmeans.bin[[1]]$a1)
                         ,R2=c( results.kmeans.bin[[1]]$R2))
      
      if(Ngroups.length>1){
        for (n in 2:Ngroups.length){
          d1<-data.frame(Method=c(paste("kmeans N=",Ngroups[n], sep="")), 
                         Bias=c(results.kmeans.bin[[n]]$bias)
                         ,MSE=c( results.kmeans.bin[[n]]$mse)
                         ,a0=c(  results.kmeans.bin[[n]]$a0)
                         ,a1=c(results.kmeans.bin[[n]]$a1)
                         ,R2=c( results.kmeans.bin[[n]]$R2))
          results=rbind(results, d1)
        }
      } 
      
    }
    
    
    dat1$pred.t1<-predicted.treat.1
    dat1$pred.t0<-predicted.treat.0
    dat1$benefit<-dat1$pred.t1-dat1$pred.t0
    glm1=glm(dat1$Y~-1+dat1$pred.t0+dat1$t:dat1$benefit, family=binomial)
    c1=format(round(coef(glm1)[1], 2), nsmall = 2)
    c2=format(round(coef(glm1)[2], 2), nsmall = 2)
    R2=format(round(rsq.n(glm1), 2), nsmall = 2)
    results.reg=data.frame("Intercept"=c1, "Slope"=c2, "Nagelkerke R2"=R2)
    rownames(results.reg)=c("Estimates") 
    
    results.all=list("kmeans"=results, "c.benefit"=results2,
                     "regression.for.benefit"=results.reg)
  }
  
  return(results.all) 
}















