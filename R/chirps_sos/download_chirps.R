SaveDir<-"/home/jovyan/common_files/chirps/raw"


DownloadCHIRPS<-function(StartYear=1980,EndYear=2021,EndDay=136,SaveDir,quiet){
  
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
    
    for(DEKAD in 1:ENDDAY){
      
      
      if(quiet){
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("Downloading file: ",DAY,"/",YEAR))
        flush.console()
      }
      
      
      DATE<-as.Date(paste0(YEAR,"-",DAY),format="%Y-%j")
      
      DAY<-format(DATE,"%d")
      MONTH<-format(DATE,"%m")
      
      FILE<-paste0("chirps-v2.0.",YEAR,".",MONTH,".",DAY,".tif.gz")
      
      URL<-paste0(URLmaster,YEAR,"/",FILE)
      destfile<-paste0(SaveDir,"/",FILE)
      if(!file.exists(destfile)){
        download.file(URL, destfile,quiet=quiet)
      }
    }
  }
}

options(timeout = 150000)
DownloadCHIRPS(StartYear=2019,EndYear = 2021,EndDay=365,SaveDir = SaveDir,quiet=T )

