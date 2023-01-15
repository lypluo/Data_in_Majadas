##################################################
###Aim:To tidy the processed PhenoCam VIs from Majadas
#-->tidy the data in order to extrat PTDs for Richard:Jan, 2023
##################################################
library(dplyr)
library(lubridate)
#-----------------
#(1)load the data from the hard drive
#-----------------
load.path<-"F:/MPI_PC/MPI-BGC/HPC/M-drive/Golden_Datasets_preparation/Data_for_use/Final_Datasets/Data_and_flag/Ecosystem/VIs/max/"
#load the North, South, Main's data:
load(paste0(load.path,"M_Eco.RDA"))
MEco_VIs<-temp;rm(temp)
load(paste0(load.path,"N_Eco.RDA"))
NEco_VIs<-temp;rm(temp)
load(paste0(load.path,"S_Eco.RDA"))
SEco_VIs<-temp;rm(temp)

##calibrate the data(Jan,2023):
#a.North tower (NT) data in 2020-->set the data after July,1 equals to data in July,1
t1<-NEco_VIs[NEco_VIs$Date>as.POSIXct("2019-01-01"),]
plot(t1$Date,t1$GCC.max_gapf)
abline(v=as.POSIXct("2020-07-01"))
#adjust the GCC data after 2020-07-01
pos<-match(as.Date("2020-07-01"),as.Date(NEco_VIs$Date))
GCC_sub<-NEco_VIs$GCC.max_gapf[pos]
NEco_VIs<-NEco_VIs %>%
  mutate(GCC.max_gapf = ifelse(Date>as.POSIXct("2020-07-01"),GCC_sub,GCC.max_gapf))

#b.Main tower (CT) data in H2016:
t2<-MEco_VIs[MEco_VIs$Date<as.POSIXct("2016-12-31"),]
plot(t2$Date,t2$GCC.max_gapf)
#in CT, the data before 2015-12 are missing-->keep the data like this

#-----------------
#(2)merge the calendar year data
#-----------------
#merge the VIs in different towers:
MEco_VIs$flag<-rep("CT",nrow(MEco_VIs))
NEco_VIs$flag<-rep("NT",nrow(NEco_VIs))
SEco_VIs$flag<-rep("NPT",nrow(SEco_VIs))
#
Eco_VIs_1<-rbind(MEco_VIs,NEco_VIs)
Eco_PhenoCam_VIs<-rbind(Eco_VIs_1,SEco_VIs)

#save the data:
save.path<-"./data/phenocam_VIs_forRichard/"
save(Eco_PhenoCam_VIs,file = paste0(save.path,"Eco_PhenoCam_VIs.RDA"))

#-----------------
#(3)also add the hydro-year data info to data:
#-----------------
#It needs to note that the GCC from CT tower is not available until Dec, 2015
# because of white balance issue..
#tidy the hydro-year info refer Richard's definition(but change the start month of the season):
#from September to next year August->Hydro2015:Sep,2014-Aug,2015:
hydro_fun<-function(df){
  # df<-MEco_VIs

  ##
  df_hydro<-df%>%
    mutate(Calyear=year(Date),
           Month=month(Date))%>%
    mutate(Hydroyear=case_when(Month<=8 ~ paste0("H",Calyear),
                               Month>8 ~ paste0("H",Calyear+1)))
  return(df_hydro)
}
#
Hydro_MEco_VIs<-hydro_fun(MEco_VIs)
Hydro_NEco_VIs<-hydro_fun(NEco_VIs)
Hydro_SEco_VIs<-hydro_fun(SEco_VIs)
#
Hydro_Eco_VIs_1<-rbind(Hydro_MEco_VIs,Hydro_NEco_VIs)
Hydro_Eco_PhenoCam_VIs<-rbind(Hydro_Eco_VIs_1,Hydro_SEco_VIs)

#save the data:
save.path<-"./data/phenocam_VIs_forRichard/"
save(Hydro_Eco_PhenoCam_VIs,file = paste0(save.path,"Hydro_Eco_PhenoCam_VIs.RDA"))

