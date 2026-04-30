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
All_ED50 <- get_all_eds_by_grouping_property(models = Allmodels)
ROUNDS=3
BIC_df=data.frame(
  Round = 0,
  Genotype=names(Allmodels),
  Nt=unlist(lapply(Allmodels,function(m){return(nrow(m$data))})),
  BIC=unlist(lapply(Allmodels,BIC)))
BIC_df=BIC_df %>% left_join(All_ED50[,c("GroupingProperty","ED50")],by=join_by("Genotype"=="GroupingProperty"))
ED50dt_df=PPpoints %>% 
  left_join(All_ED50[,c("GroupingProperty","ED50")],by=join_by("Genotype"=="GroupingProperty")) %>% 
  mutate(ED50dT=Temperature-ED50) %>% 
  group_by(Genotype) %>% 
  summarize(RMS_ED50dT=sqrt(mean(ED50dT^2)))
BIC_df=BIC_df %>% left_join(ED50dt_df,by="Genotype")

for (round in 1:ROUNDS){
  for (r in 3:7){
    #rarify subset
    this_PPp=PPpoints %>% group_by(Genotype) %>% sample_n(size=r)
    #Fit DRMS
    this_models <- fit_drms(this_PPp, grouping_properties, drm_formula, is_curveid = TRUE)
    #Get ED50 fit
    this_ED50 <- get_all_eds_by_grouping_property(models = this_models)
    #GET BIC
    this_BIC_df=data.frame(
      Round=round,
      Genotype=names(this_models),
      Nt=unlist(lapply(this_models,function(m){return(nrow(m$data))})),
      BIC=unlist(lapply(this_models,BIC)))
    this_BIC_df=this_BIC_df %>% left_join(this_ED50[,c("GroupingProperty","ED50")],by=join_by("Genotype"=="GroupingProperty"))
    this_ED50dt_df=this_PPp %>% 
      left_join(this_ED50[,c("GroupingProperty","ED50")],by=join_by("Genotype"=="GroupingProperty")) %>% 
      mutate(ED50dT=Temperature-ED50) %>% 
      group_by(Genotype) %>% 
      summarize(RMS_ED50dT=sqrt(mean(ED50dT^2)))
    this_BIC_df=this_BIC_df %>% left_join(this_ED50dt_df,by="Genotype")
    BIC_df=rbind(BIC_df,this_BIC_df)
    print(paste0("Round ",round," of ",ROUNDS,". Rarify N=",r))
  }
}
BIC_df=BIC_df %>% distinct()
BIC_df %>% ggplot(aes(x=factor(Nt),y=BIC))+geom_boxplot()
