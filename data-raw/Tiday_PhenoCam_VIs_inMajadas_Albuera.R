##################################################
###Aim:To tidy the processed PhenoCam VIs from Majadas(+Albuera)
#-->tidy the data for Pilar
##################################################

#-------------------
#1)load the data-->processed before in MPI
#-------------------
load.path<-"./data-raw/Raw_data/PhenoCam_VIs/Tree_and_Grass/Filtered_VIs/"
VI.names<-list.files(load.path)

#-------------------
#2)tidy the data for Grass
#-------------------
#load the GRASS ROIs data from different sites:
pos_GRA<-grep("_GRA",VI.names)
#
df.GRA_VIs<-c()
for (i in 1:length(pos_GRA)) {
  site.info<-VI.names[pos_GRA][i]
  if (substr(site.info,1,1)=="A"){
    sitename<-"Albuera"}
  if (substr(site.info,1,1)=="M"){
    sitename<-"Main"}
  if (substr(site.info,1,1)=="N"){
    sitename<-"North"}
  if (substr(site.info,1,1)=="S"){
    sitename<-"South"}
  #load the VIs
  load(paste0(load.path,site.info))
  temp$sitename<-rep(sitename,nrow(temp))
  #
  df.GRA_VIs<-rbind(df.GRA_VIs,temp)
  rm(temp)
}

#-------------------
#3)tidy the data for Trees
#-------------------
#load the GRASS ROIs data from different sites:
pos_EBF<-grep("_EBF",VI.names)
#
df.EBF_VIs<-c()
for (i in 1:length(pos_EBF)) {
  site.info<-VI.names[pos_EBF][i]
  if (substr(site.info,1,1)=="A"){
    sitename<-"Albuera"}
  if (substr(site.info,1,1)=="M"){
    sitename<-"Main"}
  if (substr(site.info,1,1)=="N"){
    sitename<-"North"}
  if (substr(site.info,1,1)=="S"){
    sitename<-"South"}
  #load the VIs
  load(paste0(load.path,site.info))
  temp$sitename<-rep(sitename,nrow(temp))
  #
  df.EBF_VIs<-rbind(df.EBF_VIs,temp)
  rm(temp)
}

#save the data:
save.path<-"./data/phenocam_VIs/Based_on_large_ROIs/"
save(df.GRA_VIs,file = paste0(save.path,"Grass_VI.RDA"))
save(df.EBF_VIs,file = paste0(save.path,"Tree_VI.RDA"))

#################additional######################
library(tidyverse)
library(dplyr)
library(plyr)
library(lubridate)
#-------------------
#- tidy the data for Pilar-->Majadas
#-------------------
df.GRA_VIs_Pilar<-df.GRA_VIs %>%
  filter(sitename!="Albuera")%>%
  select(Date:RCC.max_gapf,NDVI.max_gapf,NIRv.max_gapf,sitename)%>%
  mutate(Date=as.Date(Date),Year=year(Date))%>%
  filter(Year>=2017 & Year<=2019)%>%
  mutate(GCC = GCC.max_gapf,GCC.max_gapf=NULL,
         RCC = RCC.max_gapf,RCC.max_gapf=NULL,
         NDVI = NDVI.max_gapf,NDVI.max_gapf=NULL,
         NIRv = NIRv.max_gapf,NIRv.max_gapf=NULL)

df.EBF_VIs_Pilar<-df.EBF_VIs %>%
  filter(sitename!="Albuera")%>%
  select(Date:RCC.max_gapf,NDVI.max_gapf,NIRv.max_gapf,sitename)%>%
  mutate(Date=as.Date(Date),Year=year(Date))%>%
  filter(Year>=2017 & Year<=2019)%>%
  mutate(GCC = GCC.max_gapf,GCC.max_gapf=NULL,
         RCC = RCC.max_gapf,RCC.max_gapf=NULL,
         NDVI = NDVI.max_gapf,NDVI.max_gapf=NULL,
         NIRv = NIRv.max_gapf,NIRv.max_gapf=NULL)
#save the data to csv files:
save.path<-"./data/phenocam_VIs/Based_on_large_ROIs/For_Pilar/"
write.csv(df.GRA_VIs_Pilar,file = paste0(save.path,"VI_Grasses.csv"))
write.csv(df.EBF_VIs_Pilar,file = paste0(save.path,"VI_Trees.csv"))
#----------demonstration figures-----------
library(ggplot2)
df.GRA_VIs_Pilar%>%
  ggplot(aes(x=Date,y=GCC,col=sitename))+
  geom_point()

df.EBF_VIs_Pilar%>%
  ggplot(aes(x=Date,y=GCC,col=sitename))+
  geom_point()
