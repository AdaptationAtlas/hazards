
  save_dir<-"/home/jovyan/common_data/hobbins_ref_et/raw"
  
  DownloadRefET<-function(StartYear=1980,EndYear=2021,EndDekad=36,save_dir,quiet=T){
    
    
    if(!is.na(save_dir) & substr(save_dir,nchar(save_dir),nchar(save_dir))!="/"){
      save_dir<-paste0(save_dir,"/")
    }
    
    URLmaster<-"https://data.chc.ucsb.edu/products/Hobbins_RefET/ETos_p10_dekad_africa/tifs_gz/"
    

    if(!dir.exists(save_dir)){
      dir.create(save_dir,recursive=T)
    }
    
    
    for(YEAR in StartYear:EndYear){ # MinYear:MaxYear (Min = 1983 Max = Present)
      
      if(YEAR==EndYear){
        ENDDEKAD<-EndDekad
      }else{
        ENDDEKAD<-36
      }

      for(DEKAD in 1:ENDDEKAD){
        
        if(quiet){
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("Downloading file: ",DEKAD,"/",YEAR))
        flush.console()
        }
        
        
        if(DEKAD<10){
          DEKAD<-paste0("0",as.character(DEKAD))
        }
      
        FILE<-paste0("pet_",YEAR,DEKAD,".tif.gz")
      
      URL<-paste0(URLmaster,FILE)
      destfile<-paste0(save_dir,FILE)
      if(!file.exists(destfile)){
        download.file(URL, destfile,quiet=quiet)
      }
      }
    }
  }
  
  options(timeout = 150000)
  DownloadRefET(StartYear=1980,EndYear=2021,EndDekad=36,save_dir=save_dir,quiet=T)
  

  