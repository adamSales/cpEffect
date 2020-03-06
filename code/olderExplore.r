
##       `$T_Z$`=NA,
##     `$\\rho^2$`=NA
##   )

## stopifnot(max(matchTEtab[,'Estimate'])<0.25)
## stopifnot(min(matchTEtab[,'Estimate'])>0.15)

## for(tt in c(ttt['xirt'],ttt['xirt']/2))
##   matchTEtab <- add_case(matchTEtab,
##     Effect=ifelse(tt==ttt['xirt'],"[Pretest]","[Pretest/2]"),
##     Estimate=ancovat[1],
##     CI=printci(interval(tt,rrr['xirt'],b=ancovat[1], se=ancovat[2],df=ancova$df)),
##     `$\\rho^2$`=rrr['xirt'],
##     `$T_Z$`=tt)

## matchTEtab <- cbind(
##   ` `=c('\\multirow{4}{*}{Estimate}',rep('',3),'\\multirow{2}{*}{Sensitivity}',''),
##   matchTEtab
## )

## matchTEtab <- add_case(matchTEtab,` `='\\hline',.before=5)#,rep(NA,ncol(matchTEtab)-1)),matchTEtab[5:6,])

addtorow <- list(
  pos=list(0,0),
  command=c("&&&&\\multicolumn{2}{c}{Sensitivity Intervals}\\\\\n",
    paste(paste(names(matchTEtab),collapse='&'),'\\\\\n'))
)

print(xtable(matchTEtab,
  align='rrllccc',
  caption='Estimates, standard errors, and sensitivity analysis of the weighted average effect of hint usage on posttest scores, under different weighting schemes.',
  label='tab:matchResults'),
  include.rownames=FALSE,sanitize.text.function =function(x) x,add.to.row=addtorow,include.colnames=FALSE)


re1 <- ranef(mod1,condVar=TRUE)

tibble(
  slope=re1$classid2$pretestC+fixef(mod1)['pretestC'],
  V=summary(mod1)$coef['pretestC','Std. Error']+
                apply(attr(re1$classid2,'postVar'),3,function(x) x[2,2]),
  se=sqrt(V),
  ciL=slope-se,
  ciH=slope+se,
  )%>%
  arrange(slope)%>%
  mutate(ord=1:n())%>%
  ggplot(aes(ord,slope,ymin=ciL,ymax=ciH))+geom_point()+geom_errorbar(width=0)+geom_hline(yintercept=0)

mod1a <- glmer(everCP~pretestC+mpretest+(1|classid2),family=binomial,data=cpDat)

anova(mod1,mod1a)

mod2 <- update(mod1,.~.+(1|teachid2))

anova(mod1,mod2)

binnedplot(predict(mod2,type='response'),mod2@frame$everCP-predict(mod2,type='response'))


mod3 <- glmer(everCP~(1|classid2)+(1|teachid2),family=binomial,data=cpDat)

### issue with model fit: classrooms with 100% or 0% cp
cpDat <- cpDat%>%
  group_by(classid2)%>%
  mutate(pcp=mean(everCP,na.rm=TRUE),nstud=n())%>%
  ungroup()

byClass <- cpDat%>%group_by(schoolid2,teachid2,classid2,year)%>%
  summarize(pcp=pcp[1],nstud=nstud[1])

byClass%>%group_by(year)%>%summarize(p1=sum(pcp>.9),p0=sum(pcp<0.1))

byClass%>%mutate(pcpD=cut(pcp,c(-1,0.0999,0.9,2)))%>%boxplot(nstud~pcpD,data=.)

byClass%>%mutate(pcpD=cut(pcp,c(-1,0.0999,0.9,2)))%>%group_by(pcpD,year)%>%summarize(p1=mean(nstud==1),p5=mean(nstud<=5))

byClass%>%group_by(year)%>%summarize(n(),p1=sum(pcp>.9&nstud>=5),p0=sum(pcp<0.1&nstud>4))

data <- select(data,field_id,unit,section,classid2,teachid2,schoolid2,status,year)
data <- distinct(data)

ddd <- data%>%group_by(field_id,unit,section,year)%>%summarize(nrec=n(),nna=sum(is.na(status)))

ddd2 <- ddd%>%group_by(field_id)%>%mutate(anyna=any(nna==1))%>%ungroup()%>%filter(anyna)


cp <- data%>%group_by(field_id,teachid2,schoolid2,classid2,year)%>%summarize(nsec=n(),ncp=sum(status=='changed placement',na.rm=TRUE),nna=sum(is.na(status)))

boxplot(nsec~ncp,data=cp)

table(cp$ncp)/nrow(cp)
table(cp$nna)

cp <- mutate(cp,everCP=ncp>0)

boxplot(nsec~everCP,data=cp)

cpClass <- cpDat%>%group_by(classid2,teachid2,schoolid2,year)%>%summarize(nstud=n(),pcp=mean(everCP),pretest=mean(xirt,na.rm=TRUE))

qplot(pretest,pcp,color=teachid2,data=cpClass)+theme(legend.position="none")

cpTeach <- cp%>%group_by(teachid2,schoolid2,year)%>%summarize(nstud=n(),pcp=mean(everCP))
