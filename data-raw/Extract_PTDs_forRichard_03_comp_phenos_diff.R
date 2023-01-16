##################################################
###Aim:To extrat PTDs from PhenoCam GCC for Richard:Jan, 2023
##################################################
library(dplyr)
library(lubridate)
#-----------------
#(1)load the data
#-----------------
load.path<-"./data/phenocam_VIs_forRichard/extracted_phenos/"
load(file=paste0(load.path,"df.phenos_date.RDA"))
load(file=paste0(load.path,"df.phenos_doy.RDA"))
#hydro-doy
load(file=paste0(load.path,"df.phenos_hydrodoy.RDA"))

#update using the hydrodoy
#-----------------
#(2)compare the differences between treatments
#-----------------
pos_character<-match(c("Hyear","flag"),names(df.Phenos_tidy))
#NT(from H2016):
df_NT<-df.Phenos_tidy[df.Phenos_tidy$flag=="NT",-pos_character]
df_NT<-df_NT[-1,]
#NPT(from H2016)
df_NPT<-df.Phenos_tidy[df.Phenos_tidy$flag=="NPT",-pos_character]
df_NPT<-df_NPT[-1,]

#CT
df_CT<-df.Phenos_tidy[df.Phenos_tidy$flag=="CT",-pos_character]

#------comparison----------
#convert to numeric
df_NT<-apply(df_NT,2,as.numeric)
df_NPT<-apply(df_NPT,2,as.numeric)
df_CT<-apply(df_CT,2,as.numeric)

df_diff_NTCT<-df_NT - df_CT
df_diff_NPTCT<-df_NPT - df_CT
#
df_diff_NTNPT<-df_NT - df_NPT

