##################################################
###Aim:To extrat PTDs from PhenoCam GCC for Richard:Jan, 2023
##################################################
library(dplyr)
library(lubridate)

#-----------------
#(1)load the data
#-----------------
load.path<-"./data/phenocam_VIs_forRichard/"
load(file = paste0(load.path,"Hydro_Eco_PhenoCam_VIs.RDA"))


#-----------------
#(2)preparing the functions to extract the PTDs from the Majadas:
#-----------------
extract_PTDs_fun<-function(df,treatment,year){
  df<-Hydro_Eco_PhenoCam_VIs
  treatment<-"NT"
  year<-"H2015"

  #
  df.run<-df %>%
    dplyr::filter(flag==treatment & Hydroyear==year)
  #adjust some info in order to convenient for phenology extracxtion:
  df.run<-df.run %>%
    mutate(date=as.Date(Date),sitename=flag,year=Calyear,
           Date=NULL,flag=NULL,Calyear=NULL,
           doy=yday(date))

  #


}
