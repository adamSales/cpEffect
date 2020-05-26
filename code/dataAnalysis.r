library(lme4)
library(tidyverse)
library(arm)
library(RItools)
library(optmatch)
library(estimatr)
library(xtable)
library(lmerTest)
library(knitr)
library(kableExtra)
library(coefplot)
library(tikzDevice)


#dev_mode()
library(RItools)
#dev_mode()

select <- dplyr::select

`%nin%` <- Hmisc::`%nin%`

#################
### load, describe data
#################



#source('code/data.r')
load('data/data.RData')

#pload('../data/RANDstudyData/HSdata.RData')

dat <- bind_rows(
  read_csv('../data/RANDstudyData/H1_algebra_rcal_20121119_fieldid.csv'),
  read_csv('../data/RANDstudyData/H2_algebra_rcal_20121119_fieldid.csv')
)


cpDat$race6 <- dat$race[match(cpDat$field_id,dat$field_id)]

dataInfo <-
  with(cpDat,
    list(
      nstud=nrow(cpDat),
      nclass=n_distinct(classid2),
      nschool=n_distinct(schoolid2),
      states=unique(state),
      nYr1=sum(year==1),
      nYr2=sum(year==2),
      percentKept=round(nrow(cpDat)/sum(dat$treatment)*100),
      nInd=sum(race6=='AMERICAN INDIAN / ALASKAN NATIVE',na.rm=TRUE),
      nAsn=sum(race6=='ASIAN / PACIFIC ISLANDER',na.rm=TRUE),
      nOth=sum(race6=='OTHER RACE / MULTI-RACIAL',na.rm=TRUE),
      nRaceNA=sum(is.na(race6)),
      nclassTot=n_distinct(cpDat$classid2),
      nclass1=n_distinct(classid2[pcp==1]),
      nclass0=n_distinct(classid2[pcp==0]),
      pstud10=mean(pcp==1|pcp==0),
      nStud1=nrow(cpDat1)
    )
  )

### cp within classrooms



sink('tables/ncp.tex')
cpDat%>%
    mutate(
        ncp=ifelse(ncp>3,4,ncp),
        ncp=as.character(ncp),
        ncp=ifelse(ncp=='4','4+',ncp)
        )%>%
        group_by(year)%>%
        mutate(N=n())%>%
        ungroup()%>%
        group_by(year,ncp)%>%
        summarize(
            n=n(),
            `%`=round(n/N[1]*100)
            )%>%
  pivot_longer(c(n,`%`),names_to="variable",values_to="value")%>%
  pivot_wider(names_from=ncp,values_from=value)%>%
  rename(` `=variable)%>%
  kable("latex", booktabs = T, align = "c",caption='The number and percent of students in each study year of the dataset who were never reassigned, or reassigned once, twice, three times, or four or more times',
    label='ncp')%>%
  collapse_rows(columns = 1,latex_hline = "major")%>%
  add_header_above(c(" " = 2, "# Reassignments" = 5))%>%print()
sink()

#################
### propensity score match
#################


## psmod1 <- glmer(everCP~pretestC+mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+(pretestC|classid2)+(1|teachid2)+(1|schoolid2)+state,family=binomial,data=cpDat1)
## psmod0 <- update(psmod1,.~.-(pretestC|classid2)+(1|classid2))
## save(psmod1,psmod0,file='artifacts/psmod1.RData')

load('artifacts/psmod1.RData')

psAnova <- anova(psmod1,psmod0)
sdSlope <- sqrt(VarCorr(psmod1)$classid2[2,2])/sd(cpDat1$pretestC)
sdSlopep=psAnova$`Pr(>Chisq)`[2]
### coefplot
#coefplot(psmod1,intercept=FALSE,coefficients=c('race','sex','grade'))#,predictors=c('pretestC','mpretest','year'))

beta1.0 <- paste0(
  round(fixef(psmod0)['pretestC']/sd(cpDat1$pretestC),2),
  "$\\pm$",
  round(2*summary(psmod0)$coef['pretestC','Std. Error']/sd(cpDat1$pretestC),2))


psmodSumm <- summary(psmod1)

#ci1 <- confint(psmod1)

psmodSumm$coef['pretestC',1:2] <- psmodSumm$coef['pretestC',1:2]/sd(cpDat1$pretestC)
psmodSumm$coef['mpretest',1:2] <- psmodSumm$coef['mpretest',1:2]/sd(cpDat1$mpretest)

## psre <- ranef(psmod1)
## slopePos <- mean(psre$classid2$pretestC>0)
## slopeDat <- psre$classid2
## names(slopeDat) <- c('b0','b1')
## slopeDat <- add_case(slopeDat,b0=0,b1=0)
## slopeDat$b0 <- slopeDat$b0+fixef(psmod1)['(Intercept)']
## slopeDat$b1 <- slopeDat$b1+fixef(psmod1)['pretestC']

## slopeDat <- mutate(slopeDat,
##   id=1:n(),what=factor(c(rep("1 classroom",n()-1),"Average"))
##   )
## xrange=range(cpDat1$pretestC)

## yrange <- range(apply(slopeDat[,c(1,2)],1,function(x) x[1]+xrange*x[2]))


## ,
##   ymin=min(predict(psmod1,type='link')),ymax=max(predict(psmod1,type='link')))


## pdat <- cpDat1[rep(1,100),]
## pdat$pretestC <- seq(quantile(cpDat1$pretestC,.025),quantile(cpDat2$pretestC,0.975),length.out=100)
## pdat$mpretest <- mean(cpDat1$mpretest)

## pdat$yhat <- predict(psmod1,newdata=pdat,type='response')

## binnedplot(cpDat1$pretestC,cpDat1$everCP)





cpDat1$pscore1 <- predict(psmod1,type='link')

dist4 <- match_on(everCP~pscore1+xirt+strata(year),caliper=0.3,data=cpDat1)
match4 <- fullmatch(dist4,data=cpDat1)

cpDat1%>%
  mutate(
    cp=ifelse(everCP,'Reassigned','Not Reassigned'),
    excluded=ifelse(is.na(match4),'Excluded from match','Included')
  )%>%
  ggplot(aes(cp,pscore1,group=cp))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=excluded),size=0.5)+
  labs(y='Propensity Log Odds',x=NULL,color=NULL)+
  theme(legend.position='top',text=element_text(size=10))
ggsave('plots/psPlot1.pdf',width=3,height=3)

propExcluded11 <- mean(is.na(match4[cpDat1$everCP]))
propExcluded10 <- mean(is.na(match4[!cpDat1$everCP]))
propExcluded1 <- mean(is.na(match4))
nIncluded1 <- sum(!is.na(match4))
nIncluded11 <- sum(!is.na(match4)&cpDat1$everCP)
nIncluded10 <- sum(!is.na(match4)&!cpDat1$everCP)

summary(match4)

bal4 <- xBalance(everCP~xirt+pretestC+classid2+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+teachid2+schoolid2+state,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,`After Match`=~match4))

studOverall <- bal4$overall

bal4covs <- xBalance(everCP~mpretest+pretestC+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year,data=cpDat1,report=c("std.diffs","chisquare.test"),strata=list(`Before Match`=NULL,`After Match`=~match4))
covsbal <- bal4covs$overall['After.Match','p.value']

bal4class <- xBalance(everCP~classid2,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,`Matched`=~match4))

classbalSDs <- cbind(
  before=c(
    mean(bal4class$results[,'std.diff','Before.Match']<0.05,na.rm=TRUE),
    mean(bal4class$results[,'std.diff','Before.Match']<0.25,na.rm=TRUE)),
after=c(
    mean(bal4class$results[,'std.diff','Matched']<0.05,na.rm=TRUE),
    mean(bal4class$results[,'std.diff','Matched']<0.25,na.rm=TRUE))
)


classbal <- bal4class$overall

opar <- par()
pdf('plots/balance1.pdf')
plot(bal4covs,
  variable.labels=c(
     mpretest="Class mean",
      pretestC="Student",
      raceWhite='White',
      raceBlack='Black',
      raceHispanic='Latinx',
      sexF='F',
      sexM='M',
      "grade10+"="10+",
      grade9="9",
      spec_speced="Special Ed.",
      spec_gifted="Gifted",
      spec_esl="ESL",
      frl="FRL",
      frlMISTRUE="FRL Missing",
     year="Year"),
     colors=subwayPalette[1:2],
     groups=c("Pretest","Pretest","Race","Race","Race","Sex","Sex",   "Grade","Grade",NA,NA,NA,NA,NA,NA),
  strata.labels=c("--"="Before Match",match1="After Match")
  )
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
dev.off()

studCovs <- bal4covs$overall

stopifnot(studCovs['Before.Match','p.value']<0.001)


bal1Tab <- RItools:::flatten.xbalresult(bal4covs)$vartable

rownames(bal1Tab) <-
  map_chr(rownames(bal1Tab),
    ~switch(.,
      mpretest="Class Mean",
      pretestC="Class Centered",
      raceWhite='White',
      raceBlack='Black',
      raceHispanic='Latinx',
      sexF='F',
      sexM='M',
      "grade10+"="10+",
      grade9="9",
      spec_speced="Special Ed.",
      spec_gifted="Gifted",
      spec_esl="ESL",
      frl="FRL",
      frlMISTRUE="FRL Missing",
      year="Year")
 )


sink('tables/balance1.tex')
kable(bal1Tab,'latex',booktabs=TRUE,caption=paste0("Balance (standardized differences) on student level covariates before and after propensity score match. Omnibus p-values testing covariate balance are p<0.001 before matching and p=",round(studCovs['After.Match','p.value'],2)," after matching."),label='balance1',digits=2)%>%
pack_rows("Pretest",1,2)%>%pack_rows("Race/Ethnicity",3,5)%>%pack_rows("Sex",6,7)%>%pack_rows("Grade",8,9)%>%
  add_header_above(c(" " = 1, "Before Match" = 2, "After Match"=2))%>%print()
sink()

opar <- par()
pdf('plots/balancePlotStudents.pdf')
plot(bal4covs,colors=subwayPalette[1:2])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
dev.off()
do.call("par",opar)

######### go with match4

cpDat1$match <- match4

########## student level estimates
### raw
mod0 <- lm_robust(y_yirt~everCP,data=cpDat1)
### simple
(mod1 <- lm_robust(y_yirt~everCP,fixed_effects=match,data=cpDat1))
### w covariates
(mod2 <- update(mod1,.~.+splines::ns(xirt,5)+spec_speced+frlMIS+race,data=cpDat1))
### by year
(mod3.1 <- update(mod2,subset=year==1))
(mod3.2 <- update(mod2,subset=year==2))

### within classroom, plus covariates
mod4 <- lm_robust(y_yirt~everCP+splines::ns(xirt,5)+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS,fixed_effects=classid2,data=cpDat1)


## sensitivity analysis
source('code/hhh.r')

X <- as.data.frame(model.matrix(psmod1)[,-1])
X <- X[,-grep('pretest',names(X))]
X$pretest <- cpDat1$xirt
rrr <- rhos(cpDat1$y_yirt,X=X)

X2 <- X
X2$everCP <- cpDat1$everCP
ttt <- Tz(X=X2,treatment='everCP')

multipliers <- map2_dbl(rrr,ttt,~MEmult(abs(.y),.x,Inf,1.96))

stateAnova <- anova(lm(everCP~.,data=X2),lm(everCP~.-stateMI-stateKY-stateLA-stateCT-stateNJ,data=X2))

m1 <- lm(cpDat1$y_yirt~as.matrix(X))
m2 <- lm(cpDat1$y_yirt~as.matrix(select(X,-starts_with('state'))))
r1 <- summary(m1)$r.squared
r2 <- summary(m2)$r.squared
stateR <- ((1-r2)-(1-r1))/(1-r2)

effects <- list(mod0,mod1,mod2,mod3.1,mod3.2,mod4)
names(effects) <- c('mod0','mod1','mod2','mod3.1','mod3.2','mod4')

index <-map(effects,~grep('everCP',names(.$coefficients)))

sens <- map(effects,
  ~rbind(
    pre=interval(ttt['pretest'],rrr['pretest'],b=.$coefficients[grep('everCP',names(.$coefficients))],se=.$std.error[grep('everCP',names(.$coefficients))],df=Inf),
      #interval(ttt['xirt']/2,rrr['xirt']/2,b=.$coefficients,se=.$std.error,df=Inf),
      state=interval(sqrt(stateAnova$F)[2],stateR,b=.$coefficients[grep('everCP',names(.$coefficients))],se=.$std.error[grep('everCP',names(.$coefficients))],df=Inf))
)

printci <- function(ci) paste0('[',paste(sprintf('%.2f',ci),collapse=','),']')

matchTEtab <-
  map_dfr(
      names(effects),
      ~tibble(
        Model=switch(.,
          mod0='Raw',
          mod1='Matched',
          mod2='Match+Regression',
          mod3.1='Year 1',
          mod3.2='Year 2',
          mod4='Within-Class'
        ),
          Estimate=unname(effects[[.]]$coefficients[index[[.]]]),
          N=effects[[.]]$N,
        "Std. Error"=unname(effects[[.]]$std.error[index[[.]]]),
        CI=printci(with(effects[[.]],c(conf.low[index[[.]]],conf.high[index[[.]]]))),
       `[Pretest]`=printci(sens[[.]]['pre',]),
        `[State]`=printci(sens[[.]]['state',])
      )
  )

#sink('tables/effectsStud.tex')
matchTEtab <- matchTEtab%>%column_to_rownames("Model")
xtable(matchTEtab,caption='Estimates of the effect of reassignment without controlling for confounding (``Raw"), controlling for confounding with propensity score matching (``Matched"), with matching and further regression adjustment (``Match+Regression"), overall and separately for each year, and matching by classroom, with further regression adjustment (``Within-Class"). The table gives estimates, standard errors, 95\\% confidence intervals, and 95\\% sensitivity intervals assuming an unobserved confounder with properties similar to pretest scores (``[Pretest]") and to State (``[State]")',label="tab:effectsStud")%>%print(file="tables/effectsStud.tex",floating.environment="table*")
#sink()

########################
### heterogeneity
#######################

mod1re <- lmer(y_yirt ~ everCP + splines::ns(xirt, 5) + spec_speced + frlMIS +
                race+match+(everCP|classid2),data=cpDat1,REML=FALSE)
## coefTab <- summary(mod1re)$coef
## mod0re <- update(mod1re,.~.-(everCP|classid2)+(1|classid2))
## save(mod1re,mod0re,coefTab,file='artifacts/heterogeneityModel.RData')
load('artifacts/heterogeneityModel.RData')


reTest <- anova(mod1re,mod0re)
sdEffp=reTest$`Pr(>Chisq)`[2]
sdEff <- sqrt(VarCorr(mod1re)$classid2[2,2])
re1 <- ranef(mod1re,condVar=TRUE)

eff <- fixef(mod1re)['everCPTRUE']+re1$classid2[,'everCPTRUE']

effVar <- coefTab['everCPTRUE','Std. Error']^2+apply(attr(re1$classid2,'postVar'),3,function(x) x[2,2])
effSE <- sqrt(effVar)

stopifnot(coefTab['everCPTRUE','Pr(>|t|)']<0.001)

tibble(eff,effSE)%>%
  mutate(
    classid2=rownames(re1$classid2),
    Year=factor(map_int(classid2,~cpDat1$year[match(.,cpDat1$classid2)][1]))
  )%>%
  arrange(eff)%>%
  mutate(rownum=1:n())%>%
  ggplot(aes(rownum,eff,ymin=eff-effSE,ymax=eff+effSE,color=Year))+
  geom_point()+
  geom_errorbar(width=0)+
  geom_hline(yintercept=0)+
  ylab('Effect of Reassignment')+
  xlab("Classroom")+
  theme(legend.position="top")+
  annotate('text',25,.5,
    label=paste0(
      "Mean Effect=",round(fixef(mod1re)['everCPTRUE'],2),'(p<0.001)\n',
      "SD=",round(sdEff,2),' (p=',signif(reTest[['Pr(>Chisq)']][2],1),')'))
ggsave('plots/EffectByClassroom.pdf')


######## try to explain heterogeneity
cpDat1 <- cpDat1%>%
  group_by(classid2)%>%
    mutate(
      meanCP=mean(everCP,na.rm=TRUE),
      pretestVar=var(pretest,na.rm=TRUE)
      )%>%
      ungroup()
re1$classid2$se <- sqrt(apply(attr(re1$classid2,'postVar'),3,function(x) x[2,2]))

cpDat1$reTrt <- re1$classid2$everCPTRUE[cpDat1$classid2]
cpDat1$reTrtSE <- re1$classid2$se[cpDat1$classid2]

### just to look
cpDat1%>%
  group_by(classid2)%>%
  summarize_at(vars(reTrt,meanCP,pretestVar,reTrtSE),mean)%>%
  rename(`Prop. Ever Reassigned`=meanCP,`Pretest Variance`=pretestVar)%>%
  select(-classid2)%>%
  na.omit()%>%
  pivot_longer(-starts_with('reTrt'),names_to="what",values_to='x')%>%
  ggplot(aes(x,reTrt,ymin=reTrt-reTrtSE,ymax=reTrt+reTrtSE))+
  geom_point()+geom_errorbar()+geom_smooth(method='lm')+facet_wrap(~what,scale="free_x")+
  ylab('Classroom Random Slope (SE)')+xlab(NULL)
ggsave('plots/heterogeneityRegression.pdf',height=3, width=6)  

hetModcond <- update(mod1re,.~.+(meanCP+pretestVar)*everCP,data=cpDat1)
save(hetModcond,file='artifacts/hetModCond.RData')
sdEff2 <- sqrt(VarCorr(hetModcond)$classid2[2,2])
nn <- rownames(summary(hetModcond)$coef)
classHet <- summary(hetModcond)$coef[startsWith(nn,'everCPTRUE'),]

ciHet <- confint(hetModcond,rownames(classHet),method='Wald')


#####################
### classroom match
#######################
psmodClass3 <- glm(everCP~mpretest+grade9+frl+year+nstud+state,family=binomial,data=cpDatClass)

psmodClass4 <- randomForest(factor(everCP)~mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+nstud+state,data=cpDatClass)

cpDatClass$m20 <-  fullmatch(everCP~predict(psmodClass4,type='prob')[,2],data=cpDatClass)
cpDatClass$m2 <- fullmatch(psmodClass3,data=cpDatClass)

balDat <- full_join(cpDat2,tibble(classid2=cpDatClass$classid2,match1=cpDatClass$m20,match0=cpDatClass$m2))

balClass <- balanceTest(everCP~mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+nstud+cluster(classid2)+strata(match1)+strata(match0),data=balDat,report=c("std.diffs","z.scores","chisquare.test"))

balClass <- balanceTest(everCP~mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+nstud+state+schoolid2+cluster(classid2)+strata(match1),data=cpDat2,report=c("std.diffs","z.scores","chisquare.test"))

#plot(balClass2,colors=subwayPalette[1:2])

psmod2coef <- summary(psmodClass3)$coef
#for(vv in c('mpretest','grade9','frl','nstud')) psmod2coef[vv,1:2] <- psmod2coef[vv,1:2]/sd(cpDatClass[[vv]])


bal2Tab <- RItools:::flatten.xbalresult(balClass)$vartable

rownames(bal2Tab) <-
  map_chr(rownames(bal2Tab),
    ~switch(.,
      mpretest="Class Mean",
      pretestC="Class Centered",
      raceWhite='White',
      raceBlack='Black',
      raceHispanic='Latinx',
      sexF='F',
      sexM='M',
      "grade10+"="10+",
      grade9="9",
      spec_speced="Special Ed.",
      spec_gifted="Gifted",
      spec_esl="ESL",
      frl="FRL",
      frlMIS="FRL Missing",
      year="Year",
      nstud="# Students")
 )


sink('tables/balance2.tex')
kable(bal2Tab[,c(3,7,10,14)],'latex',booktabs=TRUE,caption=paste0("Balance (standardized differences) on student level covariates before and after classroom-level propensity score match. Omnibus p-values testing covariate balance are p<0.001 before matching and p=",round(studCovs['After.Match','p.value'],2)," after matching."),label='tab:balance1',digits=2)%>%
pack_rows("Pretest",1,2)%>%pack_rows("Race/Ethnicity",3,5)%>%pack_rows("Sex",6,7)%>%pack_rows("Grade",8,9)%>%
  add_header_above(c(" " = 1, "Before Match" = 2, "After Match"=2))%>%print()
sink()


source('code/balplot.r')
balPlot2 <- plotbal(balClass,colors=subwayPalette[1:2],
   variable.labels=c(
     mpretest="Pretest",
      pretestC="Student",
      raceWhite='White',
      raceBlack='Black',
      raceHispanic='Latinx',
      sexF='F',
      sexM='M',
      "grade10+"="10+",
      grade9="9",
      spec_speced="Special Ed.",
      spec_gifted="Gifted",
      spec_esl="ESL",
      frl="FRL",
      frlMIS="FRL Missing",
     year="Year"),
   strata.labels=c("--"="Before Match",match1="After Match")
)+
  geom_vline(xintercept=c(-0.25,0.25,-0.05,0.05),linetype='dotted')+
  theme(legend.position='top')+
  labs(color=NULL,shape=NULL)

ggsave('plots/balPlot2.pdf',width=3,height=4,plot=balPlot2)

mod1c <- lm_robust(y_yirt~everCP,fixed_effects=match1,clusters=cpDat2$classid2,data=cpDat2)
mod2c <- update(mod1c,.~.+race+sex+spec_gifted+spec_esl+frl+frlMIS+year)


#####################
### combined
#######################


fullCpDat <- bind_rows(cpDat2,cpDat1)

### simple
(mod1t <- lm_robust(y_yirt~everCP,fixed_effects=match,clusters=fullCpDat$classid2,data=fullCpDat))
### w covariates
(mod2t <- update(mod1t,.~.+splines::ns(xirt,5)+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+state))
### by year
(mod3.1t <- update(mod2t,subset=year==1))
(mod3.2t <- update(mod2t,subset=year==2))


save(dataInfo,eff,effSE,mod0,mod1,mod2,mod3.1,mod3.2,mod4,sens,mod1c,mod2c,mod1t,mod2t,studOverall,psmodSumm, match4,studOverall,studCovs,balClass,sdEff,sdEffp,sdSlope,sdSlopep,psmod2coef,propExcluded1,propExcluded11,propExcluded10,nIncluded1,printci,beta1.0,propExcluded11, propExcluded10, propExcluded1, nIncluded11, nIncluded10,covsbal,sdEff2,classHet,ciHet,file='artifacts/dataAnalysis.RData')
