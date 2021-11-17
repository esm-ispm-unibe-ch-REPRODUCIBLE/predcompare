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
#' @param bootstraps The number of bootstrap samples to be used for calculating confidence intervals.
#' @return A table with all estimated measures of performance.
#' @examples 
#'  # continuous outcome 
#'  dat0=simcont(500)$dat
#'  head(dat0)
#'  # Randomly shuffle the data
#'  dat<-dat0[sample(nrow(dat0)),]
#'  # Create random folds
#'  dat$folds <- cut(seq(1,nrow(dat)),breaks=10,labels=FALSE)
#'  
#'  # Obtain out-of-sample predictions
#'  dat.out.CV<-list()
#'  for (i in 1:10){
#'    dat.in.CV=dat[dat$folds!=i,]
#'    dat.out.CV[[i]]=dat[dat$folds==i,]
#'    dat1<-dat.out.CV[[i]]; dat1$t=1
#'    dat0<-dat.out.CV[[i]]; dat0$t=0
#'    m1=lm(data=dat.in.CV, y.observed~x1*t+x2*t)
#'    dat.out.CV[[i]]$predict.treat.1=predict(newdata=dat1, m1)# predictions in treatment
#'    dat.out.CV[[i]]$predict.treat.0=predict(newdata=dat0, m1)# predicions in control
#'  }
#'  
#'  dat.CV=dat.out.CV[[1]]
#'  for (i in 2:10){  dat.CV=rbind(dat.CV,dat.out.CV[[i]])}
#'  
#'  # assess model performance
#'  predcompare(repeats=20, Ngroups=c(5:10), 
#'              X=dat.CV[,c("x1", "x2","x3")], 
#'              Y=dat.CV$y.observed, 
#'              predicted.treat.1 = dat.CV$predict.treat.1,
#'              predicted.treat.0 = dat.CV$predict.treat.0,
#'              treat=dat.CV$t, type="continuous")
#'  
#'  
#'  # binary outcome 
#'  dat0=simbinary(800)$dat
#'  head(dat0)
#'  
#'  # Randomly shuffle the data
#'  dat<-dat0[sample(nrow(dat0)),]
#'  # Create random folds
#'  dat$folds <- cut(seq(1,nrow(dat)),breaks=10,labels=FALSE)
#'  
#'  dat.out.CV<-list()
#'  for (i in 1:10){
#'    dat.in.CV=dat[dat$folds!=i,]
#'    dat.out.CV[[i]]=dat[dat$folds==i,]
#'    dat1<-dat.out.CV[[i]]; dat1$t=1
#'    dat0<-dat.out.CV[[i]]; dat0$t=0
#'    glm1=glm(y.observed~(x1+x2+x3)*t, data=dat.in.CV, family = binomial(link = "logit"))
#'    dat.out.CV[[i]]$predict.treat.1=predict(newdata=dat1, glm1)# predictions in treatment
#'    dat.out.CV[[i]]$predict.treat.0=predict(newdata=dat0, glm1)# predicions in control
#'  }
#'  
#'  dat.CV=dat.out.CV[[1]]
#'  for (i in 2:10){  dat.CV=rbind(dat.CV,dat.out.CV[[i]])}
#'  
#'  
#'  predcompare(repeats=20, Ngroups=c(5:10), X=dat.CV[,c("x1", "x2","x3")], 
#'              Y=dat.CV$y.observed, 
#'              predicted.treat.1 = dat.CV$predict.treat.1,
#'              predicted.treat.0 = dat.CV$predict.treat.0,
#'              treat=dat.CV$t, type="binary", measure="RD", bootstraps = 50)
#'  
#'  predcompare(repeats=20, Ngroups=c(5:10), X=dat.CV[,c("x1", "x2","x3")], 
#'              Y=dat.CV$y.observed, 
#'              predicted.treat.1 = dat.CV$predict.treat.1,
#'              predicted.treat.0 = dat.CV$predict.treat.0,
#'              treat=dat.CV$t, type="binary", measure="logor", bootstraps = 50)
#' 
#'   
#'             
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
                     measure="RD", ### for binary outcomes only; either "RD" (risk difference) 
                     ### or "logit" is allowed
                     bootstraps=500 ### number of bootstrap samples to use in calculations
){
  library(Matching)
  library(Hmisc)
  library(MASS)
  Ngroups.length=length(Ngroups)
  
  ################ continuous data --------------
  if(type=="continuous"){
    cat(" Type of outcome: continuous","\n", 
        "Repeats: ", repeats, "\n"
    ) 
    #options(digits=2)
    
    dat1<-cbind(X,"Y"=c(Y),"t"=treat,"benefit"=predicted.treat.1-predicted.treat.0,
                "predicted.t1"=predicted.treat.1, "predicted.t0"=predicted.treat.0)
    
    
    ## calibration - group by benefit
    bias.reg<- c(); mse.reg<- c(); a0.reg<- c(); a1.reg<- c(); r2.reg<- c()
    for(group in 1:Ngroups.length){ 
      quantiles1 <- quantile(dat1$benefit, probs <-  seq(0, 1, 1/Ngroups[[group]]))
      g1 <- list(); predicted.reg <-  c();observed.reg <-  c();
      for (i in 1:Ngroups[[group]]) {g1[[i]] <- dat1[dat1$benefit >=  quantiles1[i] & 
                                                       dat1$benefit < quantiles1[i + 1], ]
      predicted.reg <-  c(predicted.reg, mean(g1[[i]]$benefit))
      observed.reg <-  c(observed.reg, mean(g1[[i]]$Y[g1[[i]]$t == 1]) - 
                           mean(g1[[i]]$Y[g1[[i]]$t == 0]))    }  
      compare=data.frame(predicted.reg, observed.reg)
      compare=compare[complete.cases(compare),]
      bias.reg<- c(bias.reg, mean((compare$predicted.reg-compare$observed.reg)))
      mse.reg<-  c(mse.reg, mean((compare$predicted.reg-compare$observed.reg)^2))
      lm1r<- summary(lm(compare$observed.reg~compare$predicted.reg))
      a0.reg<- c(a0.reg, coef(lm1r)[1,1])
      a1.reg<- c(a1.reg, coef(lm1r)[2,1])
      r2.reg<- c(r2.reg, lm1r$r.squared)
    }
    
    results.regr=data.frame(Method=c( paste("group by benefit N=",Ngroups[1], sep="")), 
                            Bias= bias.reg[1]
                            ,MSE=mse.reg[1]
                            ,a0=a0.reg[1]
                            ,a1=a1.reg[1]
                            ,R2=r2.reg[1])
    
    if(Ngroups.length>1){
      for (n in 2:Ngroups.length){
        d1<-data.frame(Method=c(paste("group by benefit N=",Ngroups[n], sep="")), 
                       Bias= bias.reg[n]
                       ,MSE=mse.reg[n]
                       ,a0=a0.reg[n]
                       ,a1=a1.reg[n]
                       ,R2=r2.reg[n])
        results.regr=rbind(results.regr, d1)
      }
    } 
    
    # Discrimination S1
    g11=dat1[dat1$benefit>0& dat1$t==1,]
    g12=dat1[dat1$benefit>0& dat1$t==0,]
    g13=dat1[dat1$benefit<0& dat1$t==1,]
    g14=dat1[dat1$benefit<0& dat1$t==0,]
    dat1$agree1=( sign(dat1$benefit) == sign(2*dat1$t-1))
    dat.disc=cbind("Y"=dat1$Y, "agree1"=dat1$agree1, X)
    s1=summary(lm(Y~., data=dat.disc))$coef[2,]
    discr.1=(paste(round(s1[1], digits=2), " [",round(s1[1]-1.96*s1[2], digits=2), "; ",  
                   round(s1[1]+1.96*s1[2], digits=2), "]", sep=""))
    
    
    # k-means clustering
    bias.m1=mse.m1=coeff.m11=coeff.m12=rsq.m1=c()
    results.kmeans<-list()
    percent1.s2<-list()
    for(n in 1:Ngroups.length){
      
      percent.s2<-c()
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
          observed.sign=(observed.groups1-observed.groups0)>0
          #estimated
          estimated.groups.m1=c(estimated.groups.m1,mean(dat1$benefit[dat1$group==i])) 
          estimated.sign=(estimated.groups.m1>0)
        } 
        # collate resuts
        dall1=(cbind(observed.sign, estimated.sign))
        dall1=dall1[complete.cases(dall1),]
        percent.s2=c(percent.s2, sum(dall1[,1]==dall1[,2] )/length(dall1[,1]))
        
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
      percent1.s2[[n]]=round(mean(percent.s2), digits=2)
      
      results.kmeans[[n]]<-data.frame(bias= round(median(bias.m1),digits=2),
                                      mse= round(median(mse.m1),digits=2),
                                      a0= round(median(coeff.m11),digits=2),
                                      a1= round(median(coeff.m12),digits=2),
                                      R2= round(median(rsq.m1),digits=2)  )        
    }
    
    # match by X ---
    bias.m1=mse.m1=m1.a0=m1.a1=m1.R2=c()
    percentX.s22<-c()
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
      percentX.s22=c(percentX.s22, mean(sign(data.compare$obs.benefit)==sign(data.compare$benefit.m1)))
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
    percentB.s2=c()
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
      percentB.s2=c(percentB.s2, mean(sign(data.compare$obs.benefit)==sign(data.compare$benefit.m1)))
      mse.m1=c(mse.m1, mean((data.compare$obs.benefit-data.compare$benefit.m1)^2))
      reg.m1=lm(data.compare$obs.benefit~data.compare$benefit.m1)
      m1.a0=c(m1.a0, coef(reg.m1)[1]); m1.a1=c(m1.a1, coef(reg.m1)[2]);  m1.R2=c(m1.R2, summary(reg.m1)$r.squared)
    }
    ## summary
    results.d=data.frame("Estimand, estimation method"=c("S1", paste("S2, k-means N=",Ngroups[1], sep="")), 
                         "Estimate"=c(discr.1, round(percent1.s2[[1]], digits=2)), check.names = F)
    results.after.matching.benefit=data.frame(bias=round(median(bias.m1),digits=2),
                                              mse=round(median(mse.m1),digits=2),
                                              a0= round(median(m1.a0),digits=2),
                                              a1= round(median(m1.a1),digits=2),
                                              R2= round(median(m1.R2),digits=2)                      )
    
    results=data.frame(Method=c( paste("kmeans N=",Ngroups[1], sep="")), 
                       Bias=c( results.kmeans[[1]]$bias)
                       ,MSE=c( results.kmeans[[1]]$mse)
                       ,a0=c( results.kmeans[[1]]$a0)
                       ,a1=c( results.kmeans[[1]]$a1)
                       ,R2=c( results.kmeans[[1]]$R2) )
    
    if(Ngroups.length>1){
      
      
      for (n in 2:Ngroups.length){
        d1<-data.frame(Method=c(paste("kmeans N=",Ngroups[n], sep="")), 
                       Bias=c(results.kmeans[[n]]$bias)
                       ,MSE=c( results.kmeans[[n]]$mse)
                       ,a0=c(  results.kmeans[[n]]$a0)
                       ,a1=c(results.kmeans[[n]]$a1)
                       ,R2=c( results.kmeans[[n]]$R2))
        results=rbind(results, d1)
        
        d2<-data.frame("Estimand, estimation method"=c(paste("S2, k-means N=",Ngroups[n], sep="")),
                       "Estimate"=as.character(c(round(percent1.s2[[n]], digits=2))), check.names = F)
        results.d=rbind(results.d, d2)
      }
    } 
    results.d=rbind(results.d, data.frame("Estimand, estimation method"="S2, match by covariates", "Estimate"=as.character(round(mean(percentX.s22), digits=2)), check.names = F))
    results.d=rbind(results.d, data.frame("Estimand, estimation method"="S2, match by benefit", "Estimate"=as.character(round(mean(percentB.s2), digits=2)), check.names = F))
    
    results=rbind(results.regr, results)
    
    results=rbind(results, data.frame(Method=c("match by covariates", "match by benefit"), 
                                      Bias=c( results.after.matching$bias ,results.after.matching.benefit$bias)
                                      ,MSE=c( results.after.matching$mse ,results.after.matching.benefit$mse)
                                      ,a0=c( results.after.matching$a0 ,results.after.matching.benefit$a0)
                                      ,a1=c( results.after.matching$a1 ,results.after.matching.benefit$a1)
                                      ,R2=c( results.after.matching$R2 ,results.after.matching.benefit$R2) ))
    
    
    
    
    
    dat1$predicted.treat.0<-predicted.treat.0
    lm1<-lm(dat1$Y~dat1$predicted.treat.0+dat1$t:dat1$benefit)
    sum.lm1<-summary(lm1)
    
    ci1=confint(lm1)
    c2=paste(round(summary(lm1)$coef[3],digits=2)," [",round( ci1[3,1],digits=2),"; " ,round(ci1[3,2], digits=2),"]", sep="")
    
    
    results.reg1=data.frame("Calibration slope"=c2)
    rownames(results.reg1)=c("Estimate") 
    
    results.all=list("outcome.type"=type,"discrimination"=results.d, "calibration"=results,
                     "regression.for.benefit"=results.reg1)
    
    
    
  }
  
  
  ################ binary data  --------------
  
  if(type=="binary"){
    library(rsq)
    cat(" Type of outcome: binary","\n", 
        "Repeats: ", repeats, "\n", 
        "Measure: ", measure, "\n", 
        "Number of bootstraps for caclulating confidence intervals: ", bootstraps,"\n" )   
    
    dat1<-cbind(X,"Y"=c(Y),"t"=treat, "benefit"=c(predicted.treat.1-predicted.treat.0))
    ### discrimination S1
    g21=dat1[dat1$benefit>0& dat1$t==1,]
    g22=dat1[dat1$benefit>0& dat1$t==0,]
    g23=dat1[dat1$benefit<0& dat1$t==1,]
    g24=dat1[dat1$benefit<0& dat1$t==0,]
    dat1$agree=1*( sign(dat1$benefit) == sign(2*dat1$t-1))
    dat.disc=cbind("Y"=dat1$Y, "agree"=dat1$agree, X)
    magree=glm(Y~., data=dat.disc, family = "binomial")
    d.test0=dat.disc; d.test0$agree=0
    p0=predict(magree, newdata=d.test0, type = "response")
    d.test1=dat.disc; d.test1$agree=1
    p1=predict(magree, newdata=d.test1, type = "response")
    S1mean=mean(p1)-mean(p0)
    bootdif <- function(dd) {
      boot=sample(nrow(dd), replace = T)
      db=dd[boot,]
      glmfit <- glm(Y~., data=db, family = "binomial")
      newdata <- db;newdata$agree=0
      p0 <- mean(predict(glmfit, newdata, type = "response"))
      newdata <- db;newdata$agree=1
      p1 <- mean(predict(glmfit, newdata, type = "response"))
      return(p1 - p0) }
    bootest <- unlist(lapply(1:bootstraps, function(x) bootdif(dat.disc)))
    S1=paste(round(S1mean, digits=3), 
             "[", round(quantile(bootest, c(0.025, 0.975))[1], digits=3), ";", 
             round(quantile(bootest, c(0.025, 0.975))[2], digits=3),"]", sep="") 
    
    ### S2 and c for benefit --------
    percent.1on1.ben=c()
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
      
      pred.ben.avg2=pred.ben.avg[obs.ben!=0]
      obs.ben2=obs.ben[obs.ben!=0]
      
      percent.1on1.ben=c(percent.1on1.ben, mean(sign(obs.ben2)==sign(pred.ben.avg2)))
      
      cindex <- rcorr.cens(pred.ben.avg, obs.ben)
      c.by.benefit <- c(c.by.benefit, cindex["C Index"][[1]])
      se.c.by.benefit<-c(se.c.by.benefit, cindex["S.D."][[1]]/2	)
      
      ## match by covariates
      percent.1on1.X=c()
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
      
      pred.ben.avg2=pred.ben.avg[obs.ben!=0]
      obs.ben2=obs.ben[obs.ben!=0]
      
      percent.1on1.X=c(percent.1on1.X, mean(sign(obs.ben2)==sign(pred.ben.avg2)))
      # Benefit c-statistic
      cindex2 <- rcorr.cens(pred.ben.avg, obs.ben)
      c.by.covariates <- c(c.by.covariates, cindex2["C Index"][[1]])
      se.c.by.covariates<-c(se.c.by.covariates, cindex2["S.D."][[1]]/2	)
      
    }
    
    mean.percent.1on1.ben=mean(percent.1on1.ben)
    mean.percent.1on1.X=mean(percent.1on1.X)
    results.c.f.b=c( 
      paste(round(c.by.covariates, digits=2), "[", round(c.by.covariates-1.96*se.c.by.covariates, digits=2), "; ", round(c.by.covariates+1.96*se.c.by.covariates, digits = 2), "]", sep=""), 
      paste(round(c.by.benefit, digits=2), "[", round(c.by.benefit-1.96*se.c.by.benefit, digits=2), "; ", round(c.by.benefit+1.96*se.c.by.benefit, digits = 2), "]", sep=""))
    
    
    
    logit<-function(x){l=log(x/(1-x)); return(l)}
    ### logor ---- 
    bias.m=mse.m=coeff.1=coeff.2=rsq.m=c()
    results.kmeans.bin<-list()
    if (measure=="logor"){
      
      
      # calibration - group by benefit
      bias.reg<- c(); mse.reg<- c(); a0.reg<- c(); a1.reg<- c(); r2.reg<- c()
      for(group in 1:Ngroups.length){ 
        quantiles1 <- quantile(dat1$benefit, probs <-  seq(0, 1, 1/Ngroups[[group]]))
        g1 <- list(); predicted.reg <-  c();observed.reg <-  c();
        for (i in 1:Ngroups[[group]]) {g1[[i]] <- dat1[dat1$benefit >=  quantiles1[i] & 
                                                         dat1$benefit < quantiles1[i + 1], ]
        predicted.reg <-  c(predicted.reg, mean(g1[[i]]$benefit))
        observed.reg <-  c(observed.reg, logit(mean(g1[[i]]$Y[g1[[i]]$t == 1])) - 
                             logit(mean(g1[[i]]$Y[g1[[i]]$t == 0])) )   }  
        compare=data.frame(predicted.reg, observed.reg)
        compare[complete.cases(compare),]
        compare <- compare[!is.infinite(rowSums(compare)),]
        bias.reg<- c(bias.reg, mean((compare$predicted.reg-compare$observed.reg)))
        mse.reg<-  c(mse.reg, mean((compare$predicted.reg-compare$observed.reg)^2))
        lm1r<- summary(lm(compare$observed.reg~compare$predicted.reg))
        a0.reg<- c(a0.reg,round( coef(lm1r)[1,1],2))
        a1.reg<- c(a1.reg,round( coef(lm1r)[2,1],2))
        r2.reg<- c(r2.reg, round(lm1r$r.squared,2))
      }
      
      results.regr=data.frame(Method=c( paste("group by benefit N=",Ngroups[1], sep="")), 
                              Bias= round(bias.reg[1],2)
                              ,MSE=round(mse.reg[1],2)
                              ,a0=round(a0.reg[1],2)
                              ,a1=round(a1.reg[1],2)
                              ,R2=round(r2.reg[1],2))
      
      if(Ngroups.length>1){
        for (n in 2:Ngroups.length){
          d1<-data.frame(Method=c(paste("group by benefit N=",Ngroups[n], sep="")), 
                         Bias=format(round(bias.reg[n],2), nsmall=2)
                         ,MSE=format(round(mse.reg[n],2), nsmall=2)
                         ,a0=a0.reg[n]
                         ,a1=a1.reg[n]
                         ,R2=r2.reg[n]
          )
          results.regr=rbind(results.regr, d1)
        }
      } 
      
      
      dat1$predicted.t1=expit(predicted.treat.1)
      dat1$predicted.t0=expit(predicted.treat.0)
      percent1.s2<-list()
      for(n in 1:Ngroups.length){
        percent.s2<-c()
        for (k in 1:repeats)
        {
          groups.kmean<-kmeans(x=dat1[,colnames(X)], centers = Ngroups[n], iter.max = 50)
          dat1$group<-groups.kmean$cluster
          ngroups<-as.vector(table(dat1$group))
          observed.groups1<-c()
          observed.groups0<-c()
          estimated.groups<-c()
          dall=NULL
          dall1=NULL
          for( i in 1:Ngroups[n]){
            # observed
            observed.groups1=c(observed.groups1, 
                               (logit(mean(dat1$Y[dat1$group==i & dat1$t==1]))))
            observed.groups0=c(observed.groups0, 
                               (logit(mean(dat1$Y[dat1$group==i & dat1$t==0]))))
            
            observed.sign=(observed.groups1-observed.groups0)>0
            # estimated
            estimated.groups=c(estimated.groups,mean(dat1$benefit[dat1$group==i])) 
            estimated.sign=(estimated.groups>0)
          } 
          # collate resuts
          dall1=(cbind(observed.sign, estimated.sign))
          dall1=dall1[complete.cases(dall1),]
          percent.s2=c(percent.s2, sum(dall1[,1]==dall1[,2] )/length(dall1[,1]))
          
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
        percent1.s2[[n]]=round(mean(percent.s2), digits=2)
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
      results=rbind(results.regr, results)
    }
    
    
    ### RD ---- 
    if (measure=="RD"){
      expit<-function(x){ exp(x)/(1 + exp(x))}
      dat1<-cbind(X,"Y"=Y,"t"=treat, "predicted.t1"=expit(predicted.treat.1),
                  "predicted.t0"=expit(predicted.treat.0))
      dat1$benefit<-(dat1$predicted.t1)-(dat1$predicted.t0)
      
      
      ## calibration - group by benefit
      bias.reg<- c(); mse.reg<- c(); a0.reg<- c(); a1.reg<- c(); r2.reg<- c()
      for(group in 1:Ngroups.length){ 
        quantiles1 <- quantile(dat1$benefit, probs <-  seq(0, 1, 1/Ngroups[[group]]))
        g1 <- list(); predicted.reg <-  c();observed.reg <-  c();
        for (i in 1:Ngroups[[group]]) {g1[[i]] <- dat1[dat1$benefit >=  quantiles1[i] & 
                                                         dat1$benefit < quantiles1[i + 1], ]
        predicted.reg <-  c(predicted.reg, mean(g1[[i]]$benefit))
        observed.reg <-  c(observed.reg, mean(g1[[i]]$Y[g1[[i]]$t == 1]) - 
                             mean(g1[[i]]$Y[g1[[i]]$t == 0]))    }  
        compare=data.frame(predicted.reg, observed.reg)
        compare=compare[complete.cases(compare),]
        bias.reg<- c(bias.reg, mean((compare$predicted.reg-compare$observed.reg)))
        mse.reg<-  c(mse.reg, mean((compare$predicted.reg-compare$observed.reg)^2))
        lm1r<- summary(lm(compare$observed.reg~compare$predicted.reg))
        a0.reg<- c(a0.reg, coef(lm1r)[1,1])
        a1.reg<- c(a1.reg, coef(lm1r)[2,1])
        r2.reg<- c(r2.reg, lm1r$r.squared)
      }
      
      results.regr=data.frame(Method=c( paste("group by benefit N=",Ngroups[1], sep="")), 
                              Bias= bias.reg[1]
                              ,MSE=mse.reg[1]
                              ,a0=a0.reg[1]
                              ,a1=a1.reg[1]
                              ,R2=r2.reg[1])
      
      if(Ngroups.length>1){
        for (n in 2:Ngroups.length){
          d1<-data.frame(Method=c(paste("group by benefit N=",Ngroups[n], sep="")), 
                         Bias= bias.reg[n]
                         ,MSE=mse.reg[n]
                         ,a0=a0.reg[n]
                         ,a1=a1.reg[n]
                         ,R2=r2.reg[n])
          results.regr=rbind(results.regr, d1)
        }
      } 
      
      
      # match via kmeans
      percent1.s2<-list()
      for(n in 1:Ngroups.length){
        percent.s2<-c()
        for (k in 1:repeats)
        {
          groups.kmean<-kmeans(x=dat1[,colnames(X)], centers = Ngroups[n], iter.max = 50)
          dat1$group<-groups.kmean$cluster
          ngroups<-as.vector(table(dat1$group))
          observed.groups1<-c()
          observed.groups0<-c()
          estimated.groups<-c()
          dall=NULL
          dall1=NULL
          for( i in 1:Ngroups[n]){
            # observed
            observed.groups1=c(observed.groups1, 
                               (mean(dat1$Y[dat1$group==i & dat1$t==1])))
            
            observed.groups0=c(observed.groups0, 
                               (mean(dat1$Y[dat1$group==i & dat1$t==0])))
            observed.sign=(observed.groups1-observed.groups0)>0
            # estimated
            estimated.groups=c(estimated.groups,mean(dat1$benefit[dat1$group==i])) 
            estimated.sign=(estimated.groups>0)
          } 
          # collate resuts
          dall1=(cbind(observed.sign, estimated.sign))
          dall1=dall1[complete.cases(dall1),]
          percent.s2=c(percent.s2, sum(dall1[,1]==dall1[,2] )/length(dall1[,1]))
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
        percent1.s2[[n]]=round(mean(percent.s2), digits=2)
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
        }      }          
      results=rbind(results.regr, results)
    }
    
    dat1$pred.t1<-predicted.treat.1
    dat1$pred.t0<-predicted.treat.0
    dat1$benefit<-dat1$pred.t1-dat1$pred.t0
    glm1=glm(dat1$Y~dat1$pred.t0+dat1$t:dat1$benefit, family=binomial)
    ci1=confint(glm1)
    slope.ben=paste(round(summary(glm1)$coef[3],digits=2)," [",round( ci1[3,1],digits=2),"; " ,round(ci1[3,2], digits=2),"]", sep="")
    
    
    results.reg=data.frame(slope.ben)
    rownames(results.reg)=c("Estimate") 
    
    
    results.dis=data.frame("Estimand, estimation method"=c("S1", paste("S2, k-means N=",Ngroups[1], sep="")), 
                           "Estimate"=c(S1, round(percent1.s2[[1]], digits=2)), check.names = F)
    if(Ngroups.length>1){
      for (n in 2:Ngroups.length){
        d1<-data.frame("Estimand, estimation method"=paste("S2, k-means N=",Ngroups[n], sep=""), 
                       "Estimate"=as.character(round(percent1.s2[[n]], digits=2)), check.names = F)
        results.dis=rbind(results.dis, d1)      }    }
    
    results.dis=rbind(results.dis,
                      data.frame("Estimand, estimation method"=c("S2, match by covariates", "S2, match by benefit",
                                                                 "c-for-benefit, match by covariates","c-for-benefit, match by benefit" ), 
                                 "Estimate"=c(as.character(round(mean.percent.1on1.X, digits=2)), as.character(round(mean.percent.1on1.ben, digits=2)), 
                                              results.c.f.b[1], results.c.f.b[2] ), check.names = F                 ))
    
    results.all=list("outcome.type"=type, "discrimination"=results.dis, "calibration.via.kmeans"=results, "regression.for.benefit"=results.reg)
    
    
  }
  
  return(results.all) 
}








