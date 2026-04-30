rm(list=ls())
#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(rstudioapi)
library(RColorBrewer)
library(CBASSED50) #Voolstra R package
fl=list.files(path = "./Data/",pattern = "_PAM.csv",full.names = T)
raw_df=NULL
for(i in 1:length(fl)){
  this_df=read.csv(fl[i])
  raw_df=rbind(raw_df,this_df)  
}
MC=as.vector(mandatory_columns())
names(raw_df)
points=raw_df %>% 
  dplyr::select(any_of(MC)) %>% 
  mutate(Colony=paste0(Species,"_s",formatC(Site,flag=0,width=2),"_c",formatC(Genotype,flag=0,width=2)),
         Genotype=Colony) %>% 
  dplyr::select(-Colony) 

#Pnt %>% group_by(Colony) %>% summarize(Nt=length(unique(Temperature))) %>% ggplot(aes(x=Nt))+geom_histogram()

#r process-and-validate-cbass-dataset
PPpoints <- preprocess_dataset(points)
validate_cbass_dataset(PPpoints)

grouping_properties <- c("Genotype")
drm_formula <- "Pam_value ~ Temperature"

Allmodels <- fit_drms(PPpoints, grouping_properties, drm_formula, is_curveid = TRUE)
#All_ED50 <- get_all_eds_by_grouping_property(models = Allmodels)
#Get ED50 Prediction and error
All_ED50.ci95=data.frame(Genotype=names(Allmodels),
                         ED50=unlist(lapply(Allmodels,function(m){
                           sm=summary(m)
                           return(sm$coefficients[3,1])
                         })),
                         ED50ci95=unlist(lapply(Allmodels,function(m){
                           sm=summary(m)
                           return(1.96*sm$coefficients[3,2])
                         })))

ROUNDS=max(sapply(X = 4:8,FUN = choose,n=8))
BIC_df=data.frame(
  Round = 0,
  Genotype=names(Allmodels),
  Nt=unlist(lapply(Allmodels,function(m){return(nrow(m$data))})),
  BIC=unlist(lapply(Allmodels,BIC)))
BIC_df=BIC_df %>% left_join(All_ED50.ci95,by=join_by("Genotype"))
ED50dt_df=PPpoints %>% 
  left_join(All_ED50.ci95[,c("Genotype","ED50")],by=join_by("Genotype")) %>% 
  mutate(ED50dT=Temperature-ED50) %>% 
  group_by(Genotype) %>% 
  summarize(RMS_ED50dT=sqrt(mean(ED50dT^2)))
BIC_df=BIC_df %>% left_join(ED50dt_df,by="Genotype")
#BIC_df$Nt[BIC_df$Nt==7]=7.5
#BIC_df$CIpval=NA
#BIC_df$BICpval=NA
minrare=4
maxrare=7
for (round in 1:ROUNDS){
  for (r in minrare:maxrare){
    #rarify subset
    this_PPp=PPpoints %>% group_by(Genotype) %>% sample_n(size=r)
    #Fit DRMS
    this_models <- fit_drms(this_PPp, grouping_properties, drm_formula, is_curveid = TRUE)
    
    # #Get ED50 fit
    # this_ED50 <- get_all_eds_by_grouping_property(models = this_models)
    
    #Get ED50 Prediction and error
    this_ED50.ci95=data.frame(Genotype=names(this_models),
                              ED50=unlist(lapply(this_models,function(m){
                                sm=summary(m)
                                return(sm$coefficients[3,1])
                              })),
                              ED50ci95=unlist(lapply(this_models,function(m){
                                sm=summary(m)
                                return(1.96*sm$coefficients[3,2])
                              })))
    
    #GET BIC
    this_BIC_df=data.frame(
      Round=round,
      Genotype=names(this_models),
      Nt=unlist(lapply(this_models,function(m){return(nrow(m$data))})),
      BIC=unlist(lapply(this_models,BIC)))
    
    #join 'em
    this_BIC_df=this_BIC_df %>% left_join(this_ED50.ci95,by=join_by("Genotype"))
    
    #get RMS_deltaT
    this_ED50dt_df=this_PPp %>% 
      left_join(All_ED50.ci95[,c("Genotype","ED50")],by=join_by("Genotype")) %>% 
      mutate(ED50dT=Temperature-ED50) %>% 
      group_by(Genotype) %>% 
      summarize(RMS_ED50dT=sqrt(mean(ED50dT^2)))
    this_BIC_df=this_BIC_df %>% left_join(this_ED50dt_df,by="Genotype")
    
    BIC_df=rbind(BIC_df,this_BIC_df)
    print(paste0("Round ",round," of ",ROUNDS,". Rarify N=",r))
  }
}
dim(BIC_df)
BIC_df=BIC_df %>% distinct()
dim(BIC_df)

for(r in (minrare+1):(maxrare+1)){
  modBIC=t.test(subset(BIC_df,Nt==(r))$BIC,subset(BIC_df,Nt==(r-1))$BIC)
  modCI=t.test(subset(BIC_df,Nt==(r))$ED50ci95,subset(BIC_df,Nt==(r-1))$ED50ci95)
  BIC_df$BICpval[BIC_df$Nt==r]=modBIC$p.value
  BIC_df$CIpval[BIC_df$Nt==r]=modCI$p.value
}

stepmet=BIC_df %>% group_by(Nt) %>% 
  summarize(meded50ci=median(ED50ci95,na.rm=T),
            medBIC=median(BIC,na.rm=T),
            pBAD=length(which(ED50ci95>=2))/length(ED50ci95))
CI_step_plot=BIC_df %>%
  filter(Nt>3) %>%
  mutate(ED50ci95=ifelse(ED50ci95>10,10,ED50ci95),
         SIG_CI=CIpval<0.05) %>% 
  ggplot(aes(x=factor(Nt),y=ED50ci95))+
  geom_violin(aes(fill=factor(SIG_CI)),quantile.linetype = c(2,1,2))+
  geom_jitter(width=.25,height=0,alpha=.1,pch=1)+
  geom_text(y=0.00001,hjust=0,nudge_x = 0.3,
            aes(label=paste0("CI95md: ",round(meded50ci,2),"\nPct Poor Fit = ",round(pBAD*100,1),"%")),data=stepmet)+
  theme_bw()+
  scale_fill_manual(name="Significant \nImprovement from\nPrevious Step",
                    values=c("TRUE"="gold","FALSE"="darkred","NA"="gray50"),drop=F)+
  xlab("Number of Temperature Steps in ED50 Fit")+
  ylab("Goodness of Fit\n ED50 95% Confidence Interval")+
  ggtitle("Goodness of Fit Improves with More Temperature Steps\n70 Rounds of Rarefaction: ED50 Confidence Interval (lower implies better fit)")+
  scale_y_log10(breaks=c(0.01,.1,1,10,20));CI_step_plot
BIC_step_plot=BIC_df %>%
  filter(Nt>3) %>%
  mutate(ED50ci95=ifelse(ED50ci95>10,10,ED50ci95),
         SIG_BIC=BICpval<0.05) %>% 
  ggplot(aes(x=factor(Nt),y=BIC))+
  geom_violin(aes(fill=factor(SIG_BIC)),quantile.linetype = c(2,1,2))+
  geom_jitter(width=.25,height=0,alpha=.1,pch=1)+
  geom_text(y=0.00001,hjust=0,nudge_x = 0.35,
            aes(label=paste0("BICmd: ",round(medBIC,2))),data=stepmet)+
  theme_bw()+
  scale_fill_manual(name="Significant \nImprovement from\nPrevious Step",
                    values=c("TRUE"="gold","FALSE"="darkred","NA"="gray50"),drop=F)+
  xlab("Number of Temperature Steps in ED50 Fit")+
  ylab("Goodness of Fit\nBayes Information Criterion")+
  ggtitle("Goodness of Fit Improves with More Temperature Steps:\n70 Round of Rarefaction: BIC (lower implies better fit)");BIC_step_plot

library(patchwork)

ABstepplot=CI_step_plot/BIC_step_plot;ABstepplot
sc=1
ggsave(plot = ABstepplot,filename = "./Plots/Rarefaction_StepPlots.jpg",width=sc*10.5,height=sc*8)

BIC_df %>% filter(Nt>3) %>% ggplot(aes(x=RMS_ED50dT,y=BIC))+
  geom_point()+xlim(c(0,6))+facet_wrap(factor(Nt)~.)
