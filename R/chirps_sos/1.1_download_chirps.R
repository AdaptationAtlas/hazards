SaveDir<-"D:/GIS Resources/Climate/CHIRPS"


DownloadCHIRPS<-function(StartYear=1980,EndYear=2022,EndDay=365,SaveDir,quiet){
  
  URLmaster<-"https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_daily/tifs/p05/"
  
  
  if(!dir.exists(SaveDir)){
    dir.create(SaveDir,recursive=T)
  }
  
  for(YEAR in StartYear:EndYear){ # MinYear:MaxYear (Min = 1983 Max = Present)
    
    if(YEAR==EndYear){
      ENDDAY<-EndDay
    }else{
      ENDDAY<-as.numeric(format(as.Date(paste0(YEAR,"-12-31")),"%j"))
    }
    
    # deal with leap years
    if(ENDDAY==365 & YEAR %in% seq(1978,2100,4)){
      ENDDAY<-366
    }
    
    for(DAY in 1:ENDDAY){
      
 
      DATE<-as.Date(paste0(YEAR,"-",DAY),format="%Y-%j")
      
      DAY<-format(DATE,"%d")
      MONTH<-format(DATE,"%m")
      
      FILE<-paste0("chirps-v2.0.",YEAR,".",MONTH,".",DAY,".tif.gz")
      FILE2<-paste0("chirps-v2.0.",YEAR,".",MONTH,".",DAY,".tif")
      
      URL<-paste0(URLmaster,YEAR,"/",FILE)
      URL2<-paste0(URLmaster,YEAR,"/",FILE2)
      destfile<-paste0(SaveDir,"/",FILE)
      destfile2<-paste0(SaveDir,"/",FILE2)
      
      if(!(file.exists(destfile)|file.exists(destfile2))){
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("Downloading file: ",FILE))
        flush.console()
        try(download.file(URL, destfile,quiet=quiet),silent=T)
      }
        
      if(!(file.exists(destfile)|file.exists(destfile2))){
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("Downloading file: ",FILE2))
        flush.console()
        try(download.file(URL2, destfile2,quiet=quiet))
      }
    }
  }
}

options(timeout = 150000)
DownloadCHIRPS(StartYear=2022,EndYear = 2022,EndDay=335,SaveDir = SaveDir,quiet=T)

