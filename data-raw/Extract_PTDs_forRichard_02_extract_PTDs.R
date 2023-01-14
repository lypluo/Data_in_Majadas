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
extract_PTDs_fun<-function(df,treatment,Hydroyear_extract){
  # df<-Hydro_Eco_PhenoCam_VIs
  # treatment<-"NT"
  # Hydroyear_extract<-"H2018"

  #
  df.run<-df %>%
    dplyr::filter(flag==treatment & Hydroyear==Hydroyear_extract)
  #adjust some info in order to convenient for phenology extracxtion:
  df.run<-df.run %>%
    mutate(date=as.Date(Date),sitename=flag,year=Calyear,
           Date=NULL,flag=NULL,Calyear=NULL,
           doy=yday(date))

  #source the extraction code:
  source(file=paste0("./R/pheno_extraction_fun_4Majadas.R"))
  df.Phenos<-SplinePheno_extraction(df.run,treatment,"GCC",FALSE,
                      as.numeric(substr(Hydroyear_extract,2,5)))
  SplinePheno_metrics_plot(df.Phenos,treatment,"GCC",FALSE,
                      as.numeric(substr(Hydroyear_extract,2,5)))
  #
  return(df.Phenos)
}

#check the PTDs extraction in every site:

#-----------for NT----------
df.Pheno.NT_2015<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2015")
df.Pheno.NT_2016<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2016")
df.Pheno.NT_2017<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2017")
df.Pheno.NT_2018<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2018") #need to check
df.Pheno.NT_2019<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2019")
df.Pheno.NT_2020<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2020") #need to figure!

#-----------for NPT----------
df.Pheno.NPT_2015<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2015")
df.Pheno.NPT_2016<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2016")
df.Pheno.NPT_2017<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2017")
df.Pheno.NPT_2018<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2018")
df.Pheno.NPT_2019<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2019")
df.Pheno.NPT_2020<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2020") #need to figure!


