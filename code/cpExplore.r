library(lme4)
library(lmerTest)
library(tidyverse)
library(arm)
library(RItools)
library(optmatch)
library(estimatr)
select <- dplyr::select

pload('data/cpPaper.RData')

cpDat <- data%>%filter(!is.na(status)& status!='final_or_incomplete')%>%
  select(field_id,year,state,section,unit,race,sex,grade,spec_speced,xirt,
    spec_gifted,spec_esl,frl,pretest,y_yirt,classid2,teachid2,schoolid2,status)%>%
  distinct()%>%
   group_by(field_id,year,state,race,sex,grade,spec_speced,xirt,
    spec_gifted,spec_esl,frl,pretest,y_yirt,classid2,teachid2,schoolid2)%>%
  summarize(nsec=n(),ncp=sum(status=='changed placement',na.rm=TRUE),nna=sum(is.na(status)))%>%
  droplevels()

cpDat$everCP <- cpDat$ncp>0

cpDat <- cpDat%>%group_by(classid2)%>%mutate(mpretest=mean(xirt,na.rm=TRUE))%>%ungroup()%>%mutate(pretestC=xirt-mpretest)

boxplot(pretestC~everCP,data=cpDat)

## mode imputation:
for(vv in c('race','sex','grade','spec_speced','spec_gifted','spec_esl','frl')){
    num <- is.numeric(cpDat[[vv]])
    cpDat[[paste0(vv,'MIS')]] <- is.na(cpDat[[vv]])
    cpDat[[vv]][is.na(cpDat[[vv]])] <- names(which.max(table(cpDat[[vv]])))
    if(num) cpDat[[vv]] <- as.numeric(cpDat[[vv]])

}

cpDat$grade <- factor(ifelse(cpDat$grade==9,'9','10+'))
levels(cpDat$race) <- list(White=c('WHITE NON-HISPANIC','ASIAN / PACIFIC ISLANDER'),Black=c('BLACK NON-HISPANIC','OTHER RACE / MULTI-RACIAL'),Hispanic=c('HISPANIC','AMERICAN INDIAN / ALASKAN NATIVE'))

cpDat$ncpCat <- factor(ifelse(cpDat$ncp>=4,'4+',cpDat$ncp))
cpDat$everCP <- cpDat$ncp>0

exirt <- stud[match(cpDat$field_id,stud$field_id),grep('Exirt2',names(stud))]
names(exirt) <- gsub('_0','_',names(exirt))
for(n in names(exirt)) cpDat[[n]] <- exirt[[n]]


cpDat <- cpDat%>%
  group_by(classid2)%>%
  mutate(pcp=mean(everCP,na.rm=TRUE),nstud=n())%>%
  ungroup()



### strategy 1: match students
cpDat1 <- filter(cpDat,pcp<1,pcp>0,!is.na(pretestC)) ### only estimate for ppl in classrooms w some but not all CP and with observed pretest

psmod1 <- glmer(everCP~pretestC+mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+(pretestC|classid2)+(1|teachid2)+(1|schoolid2)+state,family=binomial,data=cpDat1)
save(psmod1,file='artifacts/psmod1.RData')

psmod2 <- glmer(everCP~splines::ns(pretestC,5)+mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+(pretestC|classid2)+(1|teachid2)+(1|schoolid2)+state,family=binomial,data=cpDat1)

save(psmod2,file='artifacts/psmod2.RData')

psmod3 <- glmer(everCP~(splines::ns(pretestC,5)+mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS)*year+(pretestC|classid2)+(1|teachid2)+(1|schoolid2)+state,family=binomial,data=cpDat1)

save(psmod3,file='artifacts/psmod3.RData')

anova(psmod1,psmod2,psmod3)

psmod4 <- update(psmod1,.~.+pretestC:year)

### plot pretestC effect
mm <- model.matrix(psmod2)
means <- colMeans(mm)
ndat <- data.frame(pretestC=sort(cpDat1$pretestC))
for(vv in psmod1@frame%>%select(-everCP,-contains('pretestC'))%>%names()) ndat[[vv]] <- cpDat[[vv]][1]

plot(ndat$pretestC,predict(psmod2,newdata=ndat,type='link'))
points(ndat$pretestC,predict(psmod2,newdata=ndat,type='link'),pch=2)


binnedplot(predict(psmod1,type='response'),psmod1@frame$everCP-predict(psmod1,type='response'))

br <- binned.resids(predict(psmod1,type='response'),psmod1@frame$everCP-predict(psmod1,type='response'))
mean(abs(br$binned[,'ybar'])<=br$binned[,'2se'])


pscore1 <- predict(psmod1,type='link')
boxplot(pscore1~psmod1@frame$everCP)

dist1 <- match_on(everCP~pscore1,caliper=0.2,data=cpDat1)

match1 <- fullmatch(dist1,data=cpDat1)

print.xbal <- function(x) ftable(round(x$results,2), col.vars = c("strata", "stat"), row.vars = c("vars"))

bal1 <- xBalance(everCP~xirt+pretestC+classid2+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+teachid2+schoolid2+state,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,`After Match`=~match1))
bal1$overall

print(bal1covs <- xBalance(everCP~xirt+pretestC+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,`After Match`=~match1)))
bal1covs$overall

opar <- par()
plot(bal1covs,colors=subwayPalette[1:2])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

bal1class <- xBalance(everCP~classid2,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,`After Match`=~match1))

plot(bal1class,colors=subwayPalette[1:2])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)

(bal1teach <- xBalance(everCP~teachid2,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,`After Match`=~match1)))
plot(bal1teach,colors=subwayPalette[1:2])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

(bal1schoolyear <- xBalance(everCP~schoolid2*year,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,`After Match`=~match1)))
plot(bal1school,colors=subwayPalette[1:2])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

dist2 <- update(dist1,caliper=0.1)

match2 <- fullmatch(dist2,data=cpDat1)

bal2 <- xBalance(everCP~xirt+pretestC+classid2+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+teachid2+schoolid2+state,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match2))

bal2$overall

bal2covs <- xBalance(everCP~xirt+pretestC+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match2))

opar <- par()
plot(bal2covs,colors=subwayPalette[1:2])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

bal2class <- xBalance(everCP~classid2,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match2))

plot(bal2class,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)

(bal2teach <- xBalance(everCP~teachid2,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match2)))
plot(bal2teach,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

(bal2schoolyear <- xBalance(everCP~schoolid2*year,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match2)))
plot(bal2school,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

dist3 <- match_on(everCP~pscore1+xirt,caliper=0.3,data=cpDat1)
match3 <- fullmatch(dist3,data=cpDat1)

bal3 <- xBalance(everCP~xirt+pretestC+classid2+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+teachid2+schoolid2+state,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match3))

bal3$overall

(bal3covs <- xBalance(everCP~xirt+pretestC+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match3)))

opar <- par()
plot(bal3covs,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

bal3class <- xBalance(everCP~classid2,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match3))

plot(bal3class,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)

(bal3teach <- xBalance(everCP~teachid2,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match3)))
plot(bal3teach,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)

(bal3schoolyear <- xBalance(everCP~schoolid2*year,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match3)))
plot(bal3schoolyear,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)


dist4 <- match_on(everCP~pscore1+xirt+strata(year),caliper=0.3,data=cpDat1)
match4 <- fullmatch(dist4,data=cpDat1)

bal4 <- xBalance(everCP~xirt+pretestC+classid2+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+teachid2+schoolid2+state,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m1=~match1,m2=~match4))

bal4$overall

(bal4covs <- xBalance(everCP~xirt+pretestC+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year,data=cpDat1,report=c("std.diffs","z.scores","chisquare.test"),strata=list(`Before Match`=NULL,m3=~match3,m4=~match4)))

opar <- par()
plot(bal4covs,colors=subwayPalette[1:3])
abline(v=c(-0.25,-0.05,0.05,0.25),lty=2)
do.call("par",opar)


######### go with match4

cpDat1$match <- match4

#### estimate effects

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
source('~/Box Sync/rcode/hhh.r')

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

effects <- list(mod1,mod2,mod3.1,mod3.2,mod4)
names(effects) <- c('mod1','mod2','mod3.1','mod3.2','mod4')

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
        Estimate=unname(effects[[.]]$coefficients[index[[.]]]),
        "Std. Error"=unname(effects[[.]]$std.error[index[[.]]]),
        CI=printci(with(effects[[.]],c(conf.low[index[[.]]],conf.high[index[[.]]]))),
       `[Pretest]`=printci(sens[[.]]['pre',]),
        `[State]`=printci(sens[[.]]['state',])
      )
    )

### multilevel matching



### strategy 2: match classrooms
cpDat2 <- filter(cpDat,pcp==1|pcp==0,!is.na(pretestC))%>% ### only estimate for ppl in classrooms w some but not all CP and with observed pretest
  mutate(classid=as.numeric(classid2))%>%
  group_by(classid2)%>%
  mutate(nstud=n())%>%
  ungroup()


covs <- model.matrix(~pretestC+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS,data=cpDat2)
svars <- colnames(covs)[-1]
covs <- as.data.frame(covs[,-1])
covs <- cbind(covs,transmute(cpDat2,everCP=as.numeric(everCP),classid,state=as.numeric(state),year,mpretest=as.numeric(cut(mpretest,c(-Inf,quantile(unique(mpretest),seq(0.2,1,0.2)))))))

mmatch1 <- matchMulti(covs,
  treatment='everCP',
  school.id='classid',
  student.vars=svars,
  school.fb=list(c('state','year','mpretest'))
)
## dropped 23 trt and 66 ctl classrooms
#

covs <- model.matrix(~pretestC+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS,data=filter(cpDat2,nstud>5))
svars <- colnames(covs)[-1]
covs <- as.data.frame(covs[,-1])
covs <- cbind(covs,transmute(filter(cpDat2,nstud>5),everCP=as.numeric(everCP),classid,state=as.numeric(state),year,mpretest=as.numeric(cut(mpretest,c(-Inf,quantile(unique(mpretest),seq(0.25,1,0.25)))))))

mmatch2 <- matchMulti(covs,
  treatment='everCP',
  school.id='classid',
  student.vars=svars,
  school.fb=list(c('year','mpretest'))
)
## only drops one trt classroom

setdiff(mmatch1$school.match[,'TreatID'],mmatch2$school.match[,'TreatID'])
## 10 classrooms

cmatch1 <- cbind(rep(1:nrow(mmatch2$school.match),2),c(mmatch2$school.match[,1],mmatch2$school.match[,2]))

cpDat2$match1 <- cmatch1[match(cpDat2$classid,cmatch1[,2]),1]

balClass1 <- balanceTest(everCP~mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+state+cluster(classid2)+strata(match1),data=cpDat2,report=c("std.diffs","z.scores","chisquare.test"))
## um not so good. match did like nothing useful

cpDatClass <- cpDat2%>%
  group_by(classid2,schoolid2,mpretest,state,year)%>%
  mutate_if(~is.numeric(.)|is.logical(.), ~mean(.,na.rm=TRUE))%>%
  mutate(
    nstud=n(),
    black=mean(race=='Black',na.rm=TRUE),
    hisp=mean(race=='Hispanic',na.rm=TRUE),
    male=mean(sex=='M',na.rm=TRUE),
    grade9=mean(grade=='9',na.rm=TRUE)
  )%>%
  summarize_all(~.[1])%>%
  filter(nstud>5,state!='LA')

psmodClass1 <- glmer(everCP~mpretest+grade9+frl+year+nstud+(1|teachid2)+(1|schoolid2)+state,family=binomial,data=cpDatClass)
psmodClass2 <- arm::bayesglm(everCP~mpretest+grade9+frl+year+nstud+state,family=binomial,data=cpDatClass)

psmodClass3 <- glm(everCP~mpretest+grade9+frl+year+nstud+state,family=binomial,data=cpDatClass)


source('~/Box Sync/rcode/matchingFunctions.r')

scorePlot <- plotPS(psmodClass2)
par(mfrow=c(1,2))
plotPS(psmodClass2)
plotPS(psmodClass3)

par(mfrow=c(1,1))
summary(cpDatClass$m2 <- fullmatch(psmodClass3,data=cpDatClass))

plotMatch(cpDatClass$m2,psmodClass2)#,Zfuzz=scorePlot)

distC2 <- match_on(psmodClass3,caliper=0.3)
summary(cpDatClass$m3 <- fullmatch(distC2,data=cpDatClass))

plotMatch(cpDatClass$m3,psmodClass2,Zfuzz=scorePlot)

cpDat2 <- full_join(cpDat2,tibble(classid2=cpDatClass$classid2,match1=cpDatClass$m2,match2=cpDatClass$m3))

balClass2 <- balanceTest(everCP~mpretest+race+sex+grade+spec_speced+spec_gifted+spec_esl+frl+frlMIS+year+state+cluster(classid2)+strata(match1)+strata(match2),data=cpDat2,report=c("std.diffs","z.scores","chisquare.test"))

## pload('../data/RANDstudyData/HSdata.RData')

## dat$classid2 <- as.character(dat$classid2)
## cpDat2$classid2 <- as.character(cpDat2$classid2)
## cpDatClass$classid2 <- as.character(cpDatClass$classid2)

## ddd2 <- merge(cpDatClass[,c('classid2','everCP','m2')],dat,all.x=TRUE,all.y=FALSE)

## ddd <- dat%>%
##   filter(classid2%in%cpDatClass$classid2)%>%
##   mutate(
##     everCP=cpDat2$everCP[match(classid2,cpDat2$classid2)][1],
##     m2=cpDat2$m2[match(classid2,cpDat2$classid2)][1]
##   )

## print(modC <- lm_robust(
##   Y~everCP,clusters=classid2,fixed_effects=m2,data=ddd2))

## ddd2$frl <- as.numeric(as.character(ddd2$frl))
cpDat2$match <- paste0('c',as.character(cpDat2$match1))
#cpDat1$everCP <- as.numeric(cpDat1$everCP)
fullCpDat <- bind_rows(cpDat2,cpDat1)

### simple
(mod1c <- lm_robust(y_yirt~everCP,fixed_effects=match,clusters=fullCpDat$classid2,data=fullCpDat))
### w covariates
(mod2c <- update(mod1c,.~.+splines::ns(xirt,5)+spec_speced+frlMIS+race))
### by year
(mod3.1 <- update(mod2c,subset=year==1))
(mod3.2 <- update(mod2c,subset=year==2))


########################
### heterogeneity
#######################

mod1re <- lmer(y_yirt ~ everCP + splines::ns(xirt, 5) + spec_speced + frlMIS +
                 race+match+(everCP|classid2),data=cpDat1,REML=FALSE)
sdEff <- sqrt(VarCorr(mod1re)$classid2[2,2])

coefTab <- summary(mod1re)$coef

mod0re <- update(mod1re,.~.-(everCP|classid2)+(1|classid2))
reTest <- anova(mod1re,mod0re)

save(mod1re,mod0re,coefTab,file='artifacts/heterogeneityModel.RData')

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
  ggplot(aes(rownum,eff,ymin=eff-2*effSE,ymax=eff+2*effSE,color=Year))+
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

