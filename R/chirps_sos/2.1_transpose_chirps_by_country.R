require(data.table)
require(terra)
require(sf)
require(rworldxtra)
require(lubridate)
require(future)
require(future.apply)

#terraOptions(progress=0)

source("R/chirps_sos/sos_functions.R")


CHIRPS_dir="/home/jovyan/common_data/chirps_af/raw"
SaveDir<-"/home/jovyan/common_data/chirps_af/intermediate/array_countries"
CHIRPS_Dekad_Dir<-"/home/jovyan/common_data/chirps_af/intermediate/array_countries_dekad"

if(!dir.exists(CHIRPS_Dekad_Dir)){
  dir.create(CHIRPS_Dekad_Dir)
}

if(!dir.exists(SaveDir)){
  dir.create(SaveDir,recursive = T)
}

# List CHIRPS files
FILES<-data.frame(File=list.files(path = CHIRPS_dir,'.tif'))
FILES<-droplevels(FILES)

Dates<-unlist(lapply(strsplit(FILES$File,"v2.0."),"[[",2))
Dates<-unlist(lapply(strsplit(Dates,".tif"),"[[",1))
Dates<-as.Date(Dates,format = "%Y.%m.%d")

FILES$Year<-format(Dates,"%Y")
FILES$Month<-format(Dates,"%m")
FILES$Day<-format(Dates,"%d")
FILES$Date<-Dates
FILES$Yday<-format(Dates,"%j")
FILES$Code<-paste(FILES$Year,FILES$Yday,sep="-")

# Load a CHIRPS raster to get geometries
CHIRPSrast<-terra::rast(paste0(CHIRPS_dir,"/",FILES$File[1]))

# Load map of Africa
AfricaMap<-rworldmap::getMap(resolution = "high")
AfricaMap<-AfricaMap[AfricaMap$REGION=="Africa"&!is.na(AfricaMap$REGION),]
AfricaMap<-sf::st_as_sf(AfricaMap)
AfricaMap<-terra::vect(AfricaMap)
AfricaMap<-terra::project(AfricaMap,CHIRPSrast)
Countries<-sort(AfricaMap$ADMIN)
Countries<-Countries[!Countries %in% c("Cape Verde","Mauritius")]

Years<-min(FILES$Year):max(FILES$Year)

Countries<-Countries[!file.exists(paste0(SaveDir,"/",Countries,".tif"))]

# Split year stacks into countries ####
for(i in Years){
  # Display progress
  cat('\r                                                                          ')
  cat('\r',paste0("Loading data: ",i,"/",max(FILES$Year)))
  flush.console()

  Files<-paste0(SaveDir,"/",Countries,"-",i,".tif")
  
  # Load all rasters for a year
  if(!all(file.exists(Files))){
  CHIRPS<-terra::rast(paste0(CHIRPS_dir,"/",FILES$File[FILES$Year==i]))
  
  for(j in 1:length(Countries)){
    File<-Files[j]
    if(!file.exists(File)){
    # Display progress
    cat('\r                                                                      ')
    cat('\r',paste0("Processing year: ",i,"/",max(FILES$Year), " | ",Countries[j]))
    flush.console()
    
    Country<-AfricaMap[AfricaMap$ADMIN==Countries[j]]
    
    X<-terra::mask(terra::crop(CHIRPS,Country),Country)
   
    terra::writeRaster(X,file=File,overwrite=T)
    rm(X)
    gc()
    }

  }
  
  rm(CHIRPS)
  gc()
  }
}

# Combine country x year tifs ####
for(COUNTRY in  Countries){
  # Display progress
  cat('\r                                                                          ')
  cat('\r',paste0("Processing: ",COUNTRY))
  flush.console()
  
  File<-paste0(SaveDir,"/",COUNTRY,".tif")
  
  if(!file.exists(File)){
    
    Files<-paste0(SaveDir,"/",COUNTRY,"-",Years,".tif")
    
    X<-terra::rast(Files)
    terra::writeRaster(X,file=File,overwrite=T)
    
    rm(X)
    gc()
  }
  
}

# Delete year chunks
unlink(list.files(SaveDir,"-",full.names = T))

# Sum by dekad ####
source("R/chirps_sos/sos_functions.R")

for(COUNTRY in  Countries){
  # Display progress
  cat('\r                                                                          ')
  cat('\r',paste0("Processing: ",COUNTRY))
  flush.console()
  
  File<-paste0(CHIRPS_Dekad_Dir,"/",COUNTRY,".tif")
  
  if(!file.exists(File)){

    X<-terra::rast(paste0(SaveDir,"/",COUNTRY,".tif"))

    # Convert Dates into Dekads
    Dates<-names(X)
    Dates<-as.Date(gsub("chirps-v2.0.","",Dates),format=c("%Y.%m.%d"))
    Dekads<-SOS_Dekad(Dates,type="year")
    Years<-format(Dates,"%Y")
    Dekads<-paste0(Years,".",Dekads)
    
    X<-terra::tapp(X,Dekads,"sum")
    
    names(X)<-unique(Dekads)
    
    terra::writeRaster(X,file=File)
    
    rm(X)
    gc()
  }
  
}
