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
  source(file=paste0("./R/pheno_extraction_fun_4Majadas_incomplete.R"))
  #select different phenology extraction function:
  if(month(df.run$date[1])==9){
    df.Phenos<-SplinePheno_extraction(df.run,treatment,"GCC",FALSE,
                                      as.numeric(substr(Hydroyear_extract,2,5)))
    SplinePheno_metrics_plot(df.Phenos,treatment,"GCC",FALSE,
                             as.numeric(substr(Hydroyear_extract,2,5)))
  }
  if(month(df.run$date[1])>9){
    df.Phenos<-SplinePheno_extraction_incom(df.run,treatment,"GCC",FALSE,
                                      as.numeric(substr(Hydroyear_extract,2,5)))
    SplinePheno_metrics_plot_incom(df.Phenos,treatment,"GCC",FALSE,
                             as.numeric(substr(Hydroyear_extract,2,5)))
  }
  #
  return(df.Phenos)
}

#check the PTDs extraction in every site:

#-----------for NT----------
df.Pheno.NT_2015<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2015")
df.Pheno.NT_2016<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2016")
df.Pheno.NT_2017<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2017")
df.Pheno.NT_2018<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2018")
df.Pheno.NT_2019<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2019")
df.Pheno.NT_2020<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NT","H2020")
#
df.Pheno.NT_all<-list(df.Pheno.NT_2015,df.Pheno.NT_2016,df.Pheno.NT_2017,
                      df.Pheno.NT_2018,df.Pheno.NT_2019,df.Pheno.NT_2020)

#-----------for NPT----------
df.Pheno.NPT_2015<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2015")
df.Pheno.NPT_2016<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2016")
df.Pheno.NPT_2017<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2017")
df.Pheno.NPT_2018<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2018")
df.Pheno.NPT_2019<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2019")
df.Pheno.NPT_2020<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"NPT","H2020")
#
df.Pheno.NPT_all<-list(df.Pheno.NPT_2015,df.Pheno.NPT_2016,df.Pheno.NPT_2017,
                      df.Pheno.NPT_2018,df.Pheno.NPT_2019,df.Pheno.NPT_2020)

#-----------for CT----------
df.Pheno.CT_2016<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"CT","H2016")
df.Pheno.CT_2017<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"CT","H2017")
df.Pheno.CT_2018<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"CT","H2018")
df.Pheno.CT_2019<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"CT","H2019")
df.Pheno.CT_2020<-extract_PTDs_fun(Hydro_Eco_PhenoCam_VIs,"CT","H2020")
#
df.Pheno.CT_all<-list(df.Pheno.CT_2016,df.Pheno.CT_2017,
                       df.Pheno.CT_2018,df.Pheno.CT_2019,df.Pheno.CT_2020)

#-----------------
#(3)tidy the phenos:
#-----------------
tidy_pheno<-function(df,treatment){
  # df<-df.Pheno.CT_all
  # treatment<-"CT"
  #
  df.phenos<-c()
  for (i in 1:length(df)) {
    temp.phenos<-df[[i]]$Pheno_sum$pheno
    temp.phenos$Hyear<-df[[i]]$Hydro_year
    temp.phenos$flag<-treatment
    #
    df.phenos<-rbind(df.phenos,temp.phenos)
  }
  return(df.phenos)
}
#
df.Phenos_NT_tidy<-tidy_pheno(df.Pheno.NT_all,"NT")
df.Phenos_NPT_tidy<-tidy_pheno(df.Pheno.NPT_all,"NPT")
df.Phenos_CT_tidy<-tidy_pheno(df.Pheno.CT_all,"CT")
#
df.Phenos_tidy1<-rbind(df.Phenos_NT_tidy,df.Phenos_NPT_tidy)
df.Phenos_tidy<-rbind(df.Phenos_tidy1,df.Phenos_CT_tidy)

#------save the data accoring to hydrological year doy (doy=1-->Sep,1)----------
save.path<-"./data/phenocam_VIs_forRichard/extracted_phenos/"
save(df.Phenos_tidy,file=paste0(save.path,"df.phenos_hydrodoy.RDA"))

#change the hydrological year doy to Date and doy:
df.sel<-df.Phenos_tidy[,c(1:23,33)]
HydroY_to_CalY<-function(df,Hyear){
  # df<-df.sel[1,]
  # Hyear<-"H2015"
  #
  start_Date<-as.Date(paste0(as.numeric(substr(Hyear,2,5))-1,"-09-01"))
  df_CalDate<-as.Date(rep(start_Date,length(df)-1))+as.numeric(df[,-length(df)])
  df_CalDate<-t(as.data.frame(df_CalDate));df_CalDate[length(df)]<-Hyear
  df_CalDate<-as.data.frame(t(as.data.frame(df_CalDate)))
  names(df_CalDate)<-names(df)

  #
  df_Caldoy<-yday(as.Date(t(df_CalDate[,1:c(length(df)-1)])));df_Caldoy[length(df)]<-Hyear
  df_Caldoy<-as.data.frame(t(as.data.frame(df_Caldoy)))
  names(df_Caldoy)<-names(df)
  #
  df_out<-list(df_CalDate,df_Caldoy)
  return(df_out)
}
#tidying.....
df_pheno_date<-c()
df_pheno_doy<-c()

for (i in 1:nrow(df.sel)) {
  Hyear<-df.sel[i,length(df.sel)]
  df.temp<-HydroY_to_CalY(df.sel[i,],Hyear)

  #
  if(i==1){
    df_pheno_date<-df.temp[[1]]
    df_pheno_doy<-df.temp[[2]]
  }
  if(i>1){
    df_pheno_date<-rbind(df_pheno_date,df.temp[[1]])
    df_pheno_doy<-rbind(df_pheno_doy,df.temp[[2]])
  }
}
#also adding other phenological metrics:
df_pheno_date<-cbind(df_pheno_date,df.Phenos_tidy[,-c(1:23,33)])
df_pheno_doy<-cbind(df_pheno_doy,df.Phenos_tidy[,-c(1:23,33)])

#
save.path<-"./data/phenocam_VIs_forRichard/extracted_phenos/"
save(df_pheno_date,file=paste0(save.path,"df.phenos_date.RDA"))
save(df_pheno_doy,file=paste0(save.path,"df.phenos_doy.RDA"))
##save csv files for Richard:
write.csv(df_pheno_date,file=paste0(save.path,"df.phenos_date.csv"))
write.csv(df_pheno_doy,file=paste0(save.path,"df.phenos_doy.csv"))
