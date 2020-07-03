library(tidyverse)
load('data/cpPaper.RData')

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

cpDat <- cpDat%>%
  filter(!is.na(pretestC))%>%
  group_by(classid2)%>%
  mutate(pcp=mean(everCP,na.rm=TRUE),nstud=n())%>%
  ungroup()

cpDat1 <- filter(cpDat,pcp<1,pcp>0) ### only estimate for ppl in classrooms w some but not all CP and with

cpDat2 <- filter(cpDat,pcp==1|pcp==0)%>% ### only estimate for ppl in classrooms w some but not all CP and with observed pretest
  mutate(classid=as.numeric(classid2))%>%
  group_by(classid2)%>%
  mutate(nstud=n())%>%
  ungroup()


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

save(cpDat,cpDat1,cpDat2,cpDatClass,file='data/data.RData')
