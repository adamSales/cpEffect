library(arm)
library(optmatch)

between <- function(x,a,b) x>=a & x<=b
## calibrate rho and T

partialCor <- function(y,X,varb){
    ex <- extract(X,varb,TRUE)
    yres <- residuals(lm(y~.,data=ex$Xdot,na.action=na.exclude))
    Xres <- residuals(lm(ex$x~.,data=ex$Xdot,na.action=na.exclude))
    cor(yres,Xres,use='comp')
}

## takes logistic regression model and spits out
##treatment and covariates vectors
mod2dat <- function(PSmod){
 Zname <- as.character(PSmod$formula[2])
 Covnames <- names(coef(PSmod))[-1]
 Covnames <- intersect(Covnames,names(PSmod$data))
 X <- data.frame(PSmod$data[,Covnames])
 names(X)<-Covnames
 X$Z <- PSmod$data[[Zname]]
 X
}

rhos <- function(y,PSmod,X,treatment){

    if(!missing(PSmod))
        X <- model.frame(PSmod)[,-1] #mod2dat(PSmod)


    if(!missing(treatment)){
        if(is.integer(treatment)) names(X)[treatment] <- 'Z'
        if(is.character(treatment)) names(X)[which(names(X)==treatment)] <- 'Z'
    }


    rs <- NULL
    if(ncol(X)==2) rs=cor(X[names(X)!='Z'],y)
    else{
       for(i in 1:ncol(X)){
           if(names(X)[i]!='Z')
       	      rs <- c(rs,partialCor(y,X,i))
    }}
    names(rs) <- names(X)[names(X)!='Z']
    rs^2
}

newPS <- function(ex){
    psMod <- bayesglm(ex$Z~as.matrix(ex$Xdot),family=binomial(logit))
    print(summary(psMod))
    psMod$linear
}

newMatch <- function(X,mBayes,varb){
    X <- X[!is.na(mBayes),]
    ex <- extract(X,varb,FALSE)
    PS <- newPS(ex)
    match(ex$Z,PS)
}

extract <- function(X,varb,includeZ=T){
  varbCol <- ifelse(is.character(varb),which(names(X)==varb),varb)
  Xdot <- X[,-varbCol]
  x <- X[,varbCol]

  Z <- X$Z
  if(!includeZ)
    X$Z <- NULL
  list(Xdot=Xdot,Z=Z,x=x)
}


Tstat <- function(mod){
    coefs <- summary(mod)$coef
    coefs[2,1]/coefs[2,2]
}

Tz <- function(PSmod,X, treatment){
    if(!missing(PSmod)) X <- mod2dat(PSmod)

    if(!missing(treatment)){
        if(is.integer(treatment)) names(X)[treatment] <- 'Z'
        if(is.character(treatment)) names(X)[which(names(X)==treatment)] <- 'Z'
    }

    mod <- lm(Z~.,data=X)
    p <- ncol(X)
    tstats <- coef(mod)[2:p]/summary(mod)$coef[2:p,2]

    tstats
}



## compute senistivity intervals
CdT <- function(t,d)
    (1+(1+t^2)/(d-1))^(.5)

Rbound <- function(T,df,q=1.96)
    T^2/(T^2+q^2*CdT(T,df)^2)

MEmult <- function(T,R,df,q=1.96){
    if(R<Rbound(T,df,q))
        return(T*sqrt(R)+q*CdT(T,df)*(1-R)^(0.5))
    else
        return(sqrt(T^2+q^2*CdT(T,df)^2))
}

ME <- function(T,R,df,SEb,q=1.96)
    SEb*MEmult(T,R,df,q)

estMod <- function(mod, treatment='Z'){
    estRow <- summary(mod)$coef[treatment,]
    c(coef=estRow[1], se=estRow[2],df=mod$df)
}


## T is the T_W parameter from HHH, R is the rho^2_{y.w|zx} parameter
## mod is the outcome model
## treatment is the name of the treatment variable in the model
## the rest of the variables are for other uses
interval <- function(T,R,mod,treatment='Z',est,b,se,df,q=1.96,silent=TRUE){
    if(!missing(mod)){
        it=estMod(mod,treatment)
        b=it[1]
        se=it[2]
        df=it[3]
    }
    stopifnot(between(R,0,1))
    T <- abs(T)
    multiplier <- MEmult(T,R,df,q)
    if(missing(se)){
        if(!silent) print('Note: Multiplier')
        return(multiplier)
    }
    width <- multiplier*se
    if(missing(b)){
        if(!silent) print('Note: Width')
        return(width)
    }
    if(!silent) print('Interval')
    c(b-width,b+width)
}

sensitivityTable <- function(mod,X,Y,treatment,PSmod){
  T <- Tz(X = X,treatment = treatment,PSmod=PSmod)
  R <- rhos(y = Y,X = X,treatment = treatment,PSmod=PSmod)

  R <- R[order(names(T))]
  T <- T[order(names(T))]

  tab <- NULL
  for(i in 1:length(R))
    tab <- rbind(tab,
                 c(R=R[i],T=T[i],interval(T = T[i],R = R[i],mod = mod,treatment = treatment)))
  worst <- c(T=T[which.max(abs(T))],R=max(R))
  tab <- rbind(tab,
               c(R=worst['R'],T=worst[1],interval(R=worst['R'],T=worst[1],mod = mod,treatment = treatment)))
  tab <- rbind(c(R=0,T=0, interval(T = 0,R = 0,mod = mod,treatment = treatment)),tab)
  rownames(tab) <- c('No Confounding',sort(names(T)),'Worst Case')
  tab
}


