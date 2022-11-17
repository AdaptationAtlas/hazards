require(data.table)
#require(pbapply)
#require(miceadds)
require(terra)
require(sf)
install.packages("rworldxtra")
require(rworldxtra)

CHIRPS_dir<-"/home/jovyan/common_data/chirps_af/raw"
Data_dir<-"/home/jovyan/common_data/hobbins_ref_et/raw"
SaveDir<-"/home/jovyan/common_data/hobbins_ref_et/intermediate/array_countries"

if(!dir.exists(SaveDir)){
  dir.create(SaveDir,recursive = T)
}

# Unzip files
if(F){
  
  FILES<-data.table(FileName=list.files(Data_dir,".gz"))
  FILES<-FILES[!grepl("_1980",FileName)]
  
  FullPath<-list.files(Data_dir,".gz",full.names = T)
  FullPath<-FullPath[!grepl("_1980",FullPath)]
  
  FILES[,Year:=substr(FileName,5,5+3)
  ][,Dekad:=substr(FileName,9,10)
  ][,FullPath:=FullPath]
  
  for(i in 1:length(FILES$FullPath)){
    FILE<-FILES$FullPath[i]
    # Display progress
    cat('\r                                                                           ')
    cat('\r',paste0("Unzipping: ",i,"/",length(FILES$FullPath)))
    flush.console()
    
    R.utils::gunzip(FILE,remove=T,temporary=F,overwrite=T)
  }
  
}

# Split filenames into tabular form
FILES<-data.table(FileName=list.files(Data_dir,".tif"))
FILES<-FILES[!grepl("_1980",FileName)]

FullPath<-list.files(Data_dir,".tif",full.names = T)
FullPath<-FullPath[!grepl("_1980",FullPath)]

FILES[,Year:=substr(FileName,5,5+3)
      ][,Dekad:=substr(FileName,9,10)
        ][,File:=FullPath]

# Load a CHIRPS raster
CHIRPSrast<-terra::rast(list.files(CHIRPS_dir,full.names = T)[1])

# Load map of Africa
AfricaMap<-rworldmap::getMap(resolution = "high")
AfricaMap<-AfricaMap[AfricaMap$REGION=="Africa"&!is.na(AfricaMap$REGION),]
AfricaMap<-sf::st_as_sf(AfricaMap)
AfricaMap<-terra::vect(AfricaMap)
AfricaMap<-terra::project(AfricaMap,CHIRPSrast)
Countries<-sort(AfricaMap$ADMIN)
Countries<-Countries[!Countries %in% c("Cape Verde","Mauritius")]

Years<-min(FILES$Year):max(FILES$Year)


for(i in Years){
  # Display progress
  cat('\r                                                                          ')
  cat('\r',paste0("Loading data: ",i,"/",max(FILES$Year)))
  flush.console()
  
  Files<-paste0(SaveDir,"/",Countries,"-",i,".tif")
  
  # Load all rasters for a year
  if(!all(file.exists(Files))){
    STACK<-terra::rast(FILES$File[FILES$Year==i])
    
    STACK<-terra::resample(STACK,CHIRPSrast)
    
    for(j in 1:length(Countries)){
      File<-Files[j]
      if(!file.exists(File)){
        # Display progress
        cat('\r                                                                      ')
        cat('\r',paste0("Processing year: ",i,"/",max(FILES$Year), " | ",Countries[j]))
        flush.console()
        
        Country<-AfricaMap[AfricaMap$ADMIN==Countries[j]]
        
        X<-terra::mask(terra::crop(STACK,Country),Country)
        
        terra::writeRaster(X,file=File,overwrite=T)
        rm(X)
        gc()
      }
      
    }
    
    rm(STACK)
    gc()
  }
}

for(COUNTRY in  Countries){
  # Display progress
  cat('\r                                                                          ')
  cat('\r',paste0("Processing: ",COUNTRY))
  flush.console()
  
  File<-paste0(SaveDir,"/",COUNTRY,".tif")
  
  if(!file.exists(File)){
    Files<-paste0(SaveDir,"/",COUNTRY,"-",Years,".tif")
    
    X<-terra::rast(Files)
    
    writeRaster(X,file=File)
    
    rm(X)
    gc()
  }
  
}


# Delete year chunks
unlink(list.files(SaveDir,"-",full.names = T))



