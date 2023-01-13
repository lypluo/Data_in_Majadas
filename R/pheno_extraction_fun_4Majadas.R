######################################################
##Aim:functions used to extract the phenological dates
##Author:Yunpeng Luo
######################################################
#Basic steps:
#1)smoothing the time series-->spline method
#2)extract the phenological dates based on amplitude(25, 50, 75% amplitude)
#-->adding one additional EOS-->2022-02-20:EOS90 (based on 90% amplitude)
#also extract slopes at the green-up and dry-down at the same time
#--------------
#(1)function to extract the phenoloigcal dates
#--------------
# load("./data/df_GPP_Meteo_andVIs.rda")
library(phenopix)
#library(zoo)
library(plyr)
#library(hydroGOF)
library(xts)
# library(spectral.methods)
library(dplyr)
library(lubridate)
##function to norm the data
norm_data<-function(x){
  mm<-quantile(x,probs=0.95,na.rm = T)
  mn<-min(x,probs=0.05,na.rm = T)
  #update in 2024, Jan-->uisng x/c(mm-mn)-->keep the original x value info
  # norm<-(x-mn)/(mm-mn)
  norm<-(x)/(mm-mn)
  return(norm)
}

##function used to find the extremum
prod.ad1<-function(deri)
{
  Nr<-0
  len<-length(deri)-1L
  for(i in seq_len(len))
  {
    prodad<-deri[i]*deri[i+1]
    if(prodad<=0)
      Nr[i]<-i+1
  }
  Nr<-na.omit(Nr)
  Nr<-Nr[Nr!='0']
  return(Nr)
}

SplinePheno_extraction<-function(data_ts,sitename,VI_name,do_norm,year_num){
  data_ts<-df.run
  sitename<-"NT"
  VI_name<-"GCC"
  do_norm<-FALSE
  year_num<-2015

  #
  proc.ts<-data_ts %>%
    dplyr::filter(sitename==sitename)
    # dplyr::filter(year==year_num)
  proc.VIs<-proc.ts[,c("sitename","date","year","doy","Hydroyear",paste0(VI_name,".max_gapf"))]

  #-----------
  #prepare to extract the phenolgoy
  #-----------
  real_time<-proc.VIs$date

  #for selected VI name
  proc.names<-paste0(VI_name,".max_gapf")
  ts<-proc.VIs[,proc.names]

  #convert the date format to zoo
  ts<-zoo(ts,order.by = c(1:length(ts)))
  plot(ts)

  #norm the data if needed
  if(do_norm==TRUE){
    ts<-norm_data(ts)
  }

  ##gapfilling and smooth the data
  df.factor.value=0.05
  #gapfilling in order to use SplineFit afterwards
  ts<-na.fill(ts,"extend")
  fitResult<-SplineFit(ts,nrep=100,uncert= FALSE, df.factor = df.factor.value)
  ts_sm<-fitResult$fit$predicted
  rm(fitResult)

  #--------------------------
  #uncertainty estimation-->do not be used in this study
  #--------------------------
  # residuals<-as.numeric(ts)-as.numeric(ts_sm)
  # sd.res<-sd(residuals,na.rm=TRUE)
  # res2<-abs(residuals)
  # res3<-res2/max(res2)
  # sign.res<-sign(residuals)
  # predicted.df<-data.frame(matrix(ncol=100,nrow=length(ts)))
  # for(a in 1:100)
  # {
  #   noise<-runif(length(as.numeric(ts)),-sd.res,sd.res)*(as.numeric(res3)*3)
  #   sign.noise<-sign(noise)ph
  #   pos.no<-which(sign.res!=sign.noise)
  #   if(length(pos.no)!=0)
  #     noise[pos.no]<-(-noise[pos.no])
  #   fit.tmp<-SplineFit((as.numeric(ts)+as.numeric(noise)),df.factor = df.factor.value,nrep = 1)
  #   predicted.df[,a]<-as.numeric(fit.tmp$fit$predicted)
  # }
  # uncertainty<-list(predicted=predicted.df,params=NULL)
  #
  # #
  # fitResult<-list('fit'=ts_sm,'uncertainty'=uncertainty)
  # Splinefit<-fitResult
  # ###plotting
  # plot(ts,col='black',type='p')
  # for(n in 1:100)
  # {
  #   points(uncertainty$predicted[,n],col='snow2')
  # }
  # points(ts_sm,col='red')
  # }
  #evaluate Spline fitting
  # library(sirad)
  # stats_pheno<-modeval(ts_sm,ts,stat=c('RMSE','EF'))
  # #RMSE
  # SplineRMSE<-stats_pheno$RMSE
  # #R square
  # SplineRsquare<-sum((ts_sm-mean(ts))^2)/sum((ts-mean(ts))^2)
  # #Nash-Sutcliffe effciency(EF modeling effiency)
  # Spline_EF<-stats_pheno$EF
  # Splinefit_evaluation<-c('RMSE'=SplineRMSE,'R^2'=SplineRsquare,'EF'= Spline_EF)

  #------------------------------------
  #start to extract important phenological dates:timing, slope, amplitude
  #------------------------------------
  Pheno_metrics<-data.frame(matrix(ncol=1,nrow=22))
  for(a in 1:1){
    ###determine two UDs,RDs,maxs,mins from ts_sm
    deriValue<-diff(ts_sm)
    dderiValue<-diff(deriValue)
    deri_temp<-as.numeric(deriValue)
    dderi_temp<-as.numeric(dderiValue)
    DateT1<-prod.ad1(deri_temp)
    DateT2<-prod.ad1(dderi_temp)
    plot(as.numeric(ts),type = 'p')
    points(as.numeric(ts_sm),col='red')

    #-----------
    #Determine interminD
    #-----------
    baseline<- min(ts_sm, na.rm = T)
    ##define the interminD in the winter(Dec-Feb)
    ##-->extend the range 15 days before Jan and 15 days after the Feb
    extendD<-15
    intermin<-min(ts_sm[DateT1[(DateT1>122-extendD)&(DateT1<=181+extendD)]])
    interminD<-match(intermin,ts_sm)
    if(is.na(interminD)){
      intermin<-min(Ppredict[c(122+15):181])      #define between Jan-15 to end of Feb
      interminD<-match(intermin,Ppredict)
    }
    if(year==2017){
      intermin<-min(Ppredict[DateT1[(DateT1>=153-extendD)&(DateT1<=181)]]) #define between Feb to Mar
      interminD<-match(intermin,Ppredict)
    }

    abline(v=interminD,col='blue',lty=2)

    #-----------
    #start to find the transtion dates
    #-----------
    #########Determine pops
    max_line1<-max(as.numeric(ts_sm[DateT1[DateT1<interminD]]))
    max_line2<-max(as.numeric(ts_sm[DateT1[(DateT1>interminD)&(DateT1<300)]]))

    pop1<-match(max_line1,ts_sm)
    pop2<-match(max_line2,ts_sm)
    abline(v=c(pop1,pop2))

    ##Determine trs_sos&eos
    mn1<-min(ts_sm[1:180],na.rm=T)
    mn2<-min(ts_sm[181:365],na.rm = T)
    ampl1<-max_line1-mn1
    ampl2<-max_line2-mn2
    ampl_full<-max(c(max_line1,max_line2),na.rm=T)-min(c(mn1,mn2),na.rm = T)
    trs_critrion1<-mn1+0.5*ampl1
    trs_critrion2<-mn2+0.5*ampl2
    ##
    abline(h=trs_critrion1,col='blue',lty=2)
    abline(h=trs_critrion2,col='blue',lty=2)
    #
    trs_sos50_1<-30+which.min(as.numeric(abs(trs_critrion1 - ts_sm[30:pop1])))-1
    trs_eos50_1<-pop1+which.min(as.numeric(abs(trs_critrion1-ts_sm[pop1:interminD])))-1
    trs_sos50_2<-interminD+which.min(as.numeric(abs(trs_critrion2-ts_sm[interminD:pop2])))-1
    trs_eos50_2<-pop2+which.min(as.numeric(abs(trs_critrion2 - ts_sm[pop2:length(ts_sm)])))-1

    #set a crition to filter the pseduo trs50
    #if the minimum residual between trs_critrion and trs_sm bigger than
    #0.1*ampl, then set the trs50<-NA
    if(min(as.numeric(abs(trs_critrion1 - ts_sm[30:pop1])))>c(ampl1*0.05)){
      trs_sos50_1<-NA
    }
    if(min(as.numeric(abs(trs_critrion1 - ts_sm[pop1:interminD])))>c(ampl1*0.05)){
      trs_eos50_1<-NA
    }
    if(min(as.numeric(abs(trs_critrion2 - ts_sm[interminD:pop2])))>c(ampl2*0.05)){
      trs_sos50_2<-NA
    }
    if(min(as.numeric(abs(trs_critrion2 - ts_sm[pop2:length(ts_sm)])))>c(ampl2*0.05)){
      trs_eos50_2<-NA
    }

    #important to set trs_eos50_1 and  trs_sos50_2 to NA if cannot find them##
    ifelse(exists('trs_eos50_1')==FALSE,trs_eos50_1<-NA,trs_eos50_1<-trs_eos50_1)
    ifelse(exists('trs_sos50_2')==FALSE,trs_sos50_2<-NA,trs_sos50_2<-trs_sos50_2)
    abline(v=c(trs_sos50_1,trs_eos50_1,trs_sos50_2,trs_eos50_2),col='blue')
    ##-->working here-->need to work later::


    #add four more threshold-based PTDs
    trs_25_sos_value<-mn_sos+0.25*ampl_sos
    trs_75_sos_value<-mn_sos+0.75*ampl_sos
    trs_25_eos_value<-mn_eos+0.25*ampl_eos
    trs_75_eos_value<-mn_eos+0.75*ampl_eos

    trs_sos25<-15+which.min(as.numeric(abs(trs_25_sos_value - ts_sm[15:pop])))
    trs_sos75<-15+which.min(as.numeric(abs(trs_75_sos_value - ts_sm[15:pop])))
    trs_eos75<-pop+which.min(as.numeric(abs(trs_75_eos_value - ts_sm[pop:length(ts_sm)])))-1
    trs_eos25<-pop+which.min(as.numeric(abs(trs_25_eos_value - ts_sm[pop:length(ts_sm)])))-1

    #set a crition to filter the pseduo trs25,trs75
    #if the minimum residual between trs_critrion and trs_sm bigger than
    #0.1*ampl, then set the trs50<-NA
    abline(h=c(trs_25_sos_value,trs_75_sos_value,
               trs_75_eos_value,trs_25_eos_value),col='blue',lty=2)
    if(min(as.numeric(abs(trs_25_sos_value - ts_sm[15:pop])))>c(ampl_sos*0.05)){
      trs_sos25<-NA
    }
    if(min(as.numeric(abs(trs_75_sos_value - ts_sm[15:pop])))>c(ampl_sos*0.05)){
      trs_sos75<-NA
    }
    if(min(as.numeric(abs(trs_75_eos_value - ts_sm[pop:length(ts_sm)])))>c(ampl_eos*0.1)){
      trs_eos75<-NA
    }
    if(min(as.numeric(abs(trs_25_eos_value - ts_sm[pop:length(ts_sm)])))>c(ampl_eos*0.1)){
      trs_eos25<-NA
    }
    abline(v=c(trs_sos25,trs_sos75,trs_eos75,trs_eos25),col='gold')

    #add EOS 90-->add in 2022-Feb,20
    trs_90_eos_value<-mn_eos+0.9*ampl_eos
    trs_eos90<-pop+which.min(as.numeric(abs(trs_90_eos_value - ts_sm[pop:length(ts_sm)])))-1
    if(min(as.numeric(abs(trs_90_eos_value - ts_sm[pop:length(ts_sm)])))>c(ampl_eos*0.1)){
      trs_eos90<-NA
    }
    abline(v=trs_eos90,col='orange')


    ##Determine UDs and RDs
    #linear regression between [sos50-15,sos50+15], take the regression slope as the prr
    #prr-->refer Filippa et al., 2016
    library(trend)
    if(is.na(trs_sos50)){
      prr<-NA;prrD<-NA;b_prr<-NA
      UD<-NA;SD<-NA
      Gslope<-NA
    }
    if(!is.na(trs_sos50)){
      data_sub<-ts_sm[c(trs_sos50-20):c(trs_sos50+20)]
      temp_slope<-sens.slope(as.numeric(data_sub),conf.level = 0.95)
      prr<-as.numeric(temp_slope$estimates)
      prrD<-trs_sos50
      b_prr<-as.numeric(ts_sm[prrD])-prr*prrD
      #
      UD<-floor((mn_sos - b_prr)/prr)
      SD<-ceiling((max_line - b_prr)/prr)
      if(UD<0){UD<-NA;}
      if(SD>pop){SD<-NA}
      abline(v=c(UD,SD),col='green')
      abline(v=prrD,col='red')
      #add additional slope
      data_sub<-ts_sm[c(trs_sos50-20):pop]
      temp_slope<-sens.slope(as.numeric(data_sub),conf.level = 0.95)
      Gslope<-as.numeric(temp_slope$estimates)
    }

    #psr
    # linear regression between [eos50-15,eos50+15], take the regression slope as the prr
    if(is.na(trs_eos50)){
      psr<-NA;psrD<-NA;b_psr<-NA
      RD<-NA;DD<-NA
      Dslope<-NA
    }
    if(!is.na(trs_eos50)){
      if(c(trs_eos50+20)>365){
        psrD<-trs_eos50
        psr<-NA;b_psr<-NA
        RD<-NA;DD<-NA
        Dslope<-NA
      }
      if(c(trs_eos50+20)<=365){
        data_sub<-ts_sm[c(trs_eos50-20):c(trs_eos50+20)]
        temp_slope<-sens.slope(as.numeric(data_sub),conf.level = 0.95)
        psr<-as.numeric(temp_slope$estimates)
        psrD<-trs_eos50
        b_psr<-as.numeric(ts_sm[psrD])-psr*psrD
        #
        RD <- floor((mn_eos - b_psr)/psr)
        DD<-ceiling((max_line-b_psr)/psr)
        abline(v=c(RD,DD),col='green')
        abline(v=psrD,col='red')
        #add additional slope
        # if(RD<=length(ts_sm)){
          ##add additional slope
          data_sub<-ts_sm[pop:c(trs_eos50+20)]
          temp_slope<-sens.slope(as.numeric(data_sub),conf.level = 0.95)
          Dslope<-as.numeric(temp_slope$estimates)
        # }
        if(RD>length(ts_sm)){
          ##add additional slope
          RD<-NA
        }

      }

    }
    #put the phenometrics together
    Pheno_metrics[,a]<- c(UD,prrD,SD,pop,DD,psrD,RD,
                          trs_sos25,trs_sos50,trs_sos75,trs_eos90,trs_eos75,trs_eos50,trs_eos25,
                          max_line,baseline,ampl_sos,ampl_eos,prr,psr,Gslope,Dslope)
  }
  Pheno_metrics<-as.data.frame(t(Pheno_metrics))
  names(Pheno_metrics)<-c("UD","prrD","SD","pop","DD","psrD","RD",
                              'trs_sos25','trs_sos50','trs_sos75',"trs_eos90",'trs_eos75','trs_eos50','trs_eos25',
                              "maxline","baseline",'ampl_sos','ampl_eos',"prr","psr","Gslope","Dslope")
  row.names(Pheno_metrics)<-site.name
  #
  #some adjust: in case the time series is not start from doy=1:
  min.doy<-min(proc.VIs$doy,na.rm=T)
  #some adjust: in case the time series is not start from doy=1:
  Pheno_metrics[1:14]<-Pheno_metrics[1:14]+min.doy-1
  #----------------------
  # metrics summary
  #---------------------------
  metrics<- c(UD,prrD,SD,pop,DD,psrD,RD)
  metrics_part1<-c(UD,pop)
  metrics_part2<-c(prrD,psrD)
  metrics_part3<-c(SD,DD)
  metrics_trs<-c(trs_sos25,trs_sos50,trs_sos75,trs_eos90,trs_eos75,trs_eos50,trs_eos25)
  names(metrics) <- c("UD","prrD","SD","pop","DD","psrD","RD")
  names(metrics_part1)<-c("UD","pop")
  names(metrics_part2)<-c("prrD","psrD")
  names(metrics_part3)<-c("SD","DD")
  names(metrics_trs)<-c('trs_sos25','trs_sos50','trs_sos75',"trs_eos90",'trs_eos75','trs_eos50','trs_eos25')
  allparas<-c("UD"=UD,"prrD"=prrD,"SD"=SD,"pop"=pop,"DD"=DD,"psrD"=psrD,"RD"=RD,
              "trs_sos25"=trs_sos25,"trs_sos50"=trs_sos50,"trs_sos75"=trs_sos75,
              "trs_eos25"=trs_eos25,"trs_eos50"=trs_eos50,"trs_eos75"=trs_eos75,"trs_eos90"=trs_eos90,
              "maxline"=max_line,"baseline"=baseline,"ampl_sos"=ampl_sos,"ampl_eos"=ampl_eos,
              "prr"= prr,"psr"=psr,"Gslope"=Gslope,"Dslope"=Dslope)

  #for data synthesis
  Pheno_sum<-c()
  Pheno_sum$allparas<-allparas
  # Pheno_sum$Splinefit<-Splinefit
  Pheno_sum$pheno<-Pheno_metrics

  #for preparing the plots
  all_ts<-data.frame(time=real_time,ts_ori=ts[,VI_name],ts_fit=ts_sm,ts=ts)
  plot_paras<-list(b_prr=b_prr,prr=prr,b_psr=b_psr,psr=psr,SD=SD,DD=DD,SD=SD,DD=DD,baseline=baseline,
                   max_line=max_line,
                   trs_critrion_sos=trs_critrion_sos,trs_critrion_eos=trs_critrion_eos)
  plot_metrics<-list(metrics=metrics,metrics_part1=metrics_part1,metrics_part2=metrics_part2,metrics_part3=metrics_part3,
                     metrics_trs=metrics_trs,plot_paras=plot_paras)
  Pheno_result<-list(Pheno_sum=Pheno_sum,all_ts=all_ts,plot_metrics=plot_metrics)
  return(Pheno_result)
}

#--------------
#(2)function to plotting the extract the results
#--------------
##Pheno_metrics visualization
# data_ts<-df.ts
# site.name<-"LS"
# VI_name<-"EVI"
# do_norm<-FALSE
# year_num<-2000
# Pheno_results<-SplinePheno_extraction(data_ts,site.name,VI_name,FALSE,year)
SplinePheno_metrics_plot<-function(Pheno_results,site.name,VI_name,do_norm,Year){
  # Pheno_results<-Pheno_result
  # site.name<-"ROS2"
  # VI_name<-"GPP"
  # Year<-2015
  # do_norm<-FALSE
  #
  #
  ts<-Pheno_results$all_ts$ts_ori
  ts_sm<-Pheno_results$all_ts$ts_fit
  deriValue<-diff(as.numeric(ts_sm))
  #add in Jan,2023
  doy=yday(Pheno_results$all_ts$time)
  ts<-zoo(ts,order.by = doy)
  ts_sm<-zoo(ts_sm,order.by = doy)
  #
  yname<-ifelse(do_norm==TRUE,"Norm-","")

  #(1)plot1:  gcc time series and its first derivitives
  # par(mfrow=c(2,1))
  # plot(ts,xlab=paste0(site.name,'_Date'),ylab=paste0(VI_name),
  #     main=paste0(site.name,' ',VI_name,' variation(spline)'),type='p',xaxt="n")
  # axis(1,at=c(0,30,61,91,122,153,181,212,242,273,303,334,365),labels=c("Sep",
  # "Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",'Sep'))
  # points(ts_sm,col='red')
  # plot(deriValue,xlab='Date',ylab=paste0(VI_name,'_deri'),col='green',xaxt="n")
  # axis(1,at=c(0,30,61,91,122,153,181,212,242,273,303,334,365),labels=c("Sep",
  # "Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",'Sep'))

  #(2)plot 2: demo plot for the extraction results
  ###########
  plot_metrics<-Pheno_results$plot_metrics
  ##adding the values:
  plot_metrics$metrics<-plot_metrics$metrics+doy[1]-1
  plot_metrics$metrics_part1<-plot_metrics$metrics_part1+doy[1]-1
  plot_metrics$metrics_part2<-plot_metrics$metrics_part2+doy[1]-1
  plot_metrics$metrics_part3<-plot_metrics$metrics_part3+doy[1]-1
  plot_metrics$metrics_trs<-plot_metrics$metrics_trs+doy[1]-1
  #
  b_prr<-plot_metrics$plot_paras$b_prr+doy[1]-1;prr<-plot_metrics$plot_paras$prr
  b_psr<-plot_metrics$plot_paras$b_psr+doy[1]-1;psr<-plot_metrics$plot_paras$psr
  SD<-plot_metrics$plot_paras$SD+doy[1]-1;DD<-plot_metrics$plot_paras$DD1+doy[1]-1
  baseline<-plot_metrics$plot_paras$baseline;max_line<-plot_metrics$plot_paras$max_line
  trs_critrion_sos<-plot_metrics$plot_paras$trs_critrion_sos+doy[1]-1
  trs_critrion_eos<-plot_metrics$plot_paras$trs_critrion_eos+doy[1]-1

  par(mfrow=c(2,1))
  ##plot1
  plot(ts,xlab=paste0(site.name,'_Spline_',Year),ylab=paste0(yname,VI_name),main=paste('phenophase-estimate'),xaxt='n',type = 'p')
  points(ts_sm,col='red')
  axis(1,at=c(0,30,61,91,122,153,181,212,242,273,303,334,365),labels=c("Jan",
  "Feb","Mar","Apr","May","Jun","Jul","Aug",'Sep',"Oct","Nov","Dec","Jan"))
  #colors <- palette()[2:5]
  colors<-c('lightgreen','blue','cyan','red','cyan','blue','orange')
  colors1<-c('lightgreen','red')
  colors2<-c('blue','blue')
  colors3<-c('cyan','cyan')
  ylons1 <- c(min(abs(ts_sm), na.rm = TRUE) * 1.01)
  ylons2 <- c(min(abs(ts_sm), na.rm = TRUE) * 1.05)
  ylons3 <- c(min(abs(ts_sm), na.rm = TRUE) * 1.1)
  abline(v = plot_metrics$metrics, col = colors)
  text(plot_metrics$metrics_part1, y = ylons1,labels = names(plot_metrics$metrics_part1),col = colors1,cex=0.75)
  text(plot_metrics$metrics_part2, y = ylons2,labels = names(plot_metrics$metrics_part2),col = colors2,cex=0.75)
  text(plot_metrics$metrics_part3, y = ylons3,labels = names(plot_metrics$metrics_part3),col = colors3,cex=0.75)
  if(!is.na(b_prr)){
    abline(a=b_prr,b=prr,col='blue',lty=2,lwd='2.5')
  }
  if(!is.na(b_psr)){
    abline(a=b_psr,b=psr,col='blue',lty=2,lwd='2.5')
  }
  abline(h=c(baseline,max_line),col='green',lty= 2,lwd='2.5')

  ##plot2
  plot(ts,xlab=paste0(site.name,'_Spline_',Year),ylab=paste0(yname,VI_name),
       main=paste('phenophase-estimate'),xaxt='n',type = 'p')
  points(ts_sm,col='red')
  axis(1,at=c(0,30,61,91,122,153,181,212,242,273,303,334,365),labels=c("Jan",
     "Feb","Mar","Apr","May","Jun","Jul","Aug",'Sep',"Oct","Nov","Dec","Jan"))
  colors<-c('lightgreen','blue','cyan','red','cyan','blue','orange')
  colors1<-c('blue','blue','blue','blue','blue','blue')
  ylons <- c(min(ts_sm, na.rm = T) * 1.01)
  ylons1 <- c(min(ts_sm, na.rm = T) *2)
  abline(v = plot_metrics$metrics, col = colors)
  text(plot_metrics$metrics, y = ylons, labels = names(plot_metrics$metrics),col = colors,cex=0.75)
  abline(v=plot_metrics$metrics_trs,col=colors1)
  text(plot_metrics$metrics_trs, y = ylons1, labels = substr(names(plot_metrics$metrics_trs),5,9),col = colors1,cex=0.75)
}

#plotting:
# SplinePheno_metrics_plot(Pheno_results,site.name,"EVI",year)