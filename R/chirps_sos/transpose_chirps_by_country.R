require(data.table)
require(pbapply)
require(R.utils)
require(miceadds)
require(terra)
require(sf)

rgdal::set_proj_search_paths("/opt/conda/share/proj")


CHIRPS_dir="/home/jovyan/common_data/chirps_af/raw"
SaveDir<-"/home/jovyan/common_data/chirps_af/intermediate/array_countries"

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
CHIRPSrast<-terra::classify(CHIRPSrast, c(-9999,NA))
terra::plot(CHIRPSrast)

# Load map of Africa
AfricaMap<-rworldmap::getMap(resolution = "high")
AfricaMap<-AfricaMap[AfricaMap$REGION=="Africa"&!is.na(AfricaMap$REGION),]
AfricaMap<-sf::st_as_sf(AfricaMap)
AfricaMap<-terra::vect(AfricaMap)
AfricaMap<-terra::project(AfricaMap,CHIRPSrast)
AfricaMap<-terra::rasterize(AfricaMap,CHIRPSrast,field="ADMIN")

# Convert to points
points<-terra::as.points(AfricaMap,values=T,na.rm=F)

# Make an index in case we need to use it to cross reference list items back to the raster
# Choose timestep in middle of sequence where missing values unlikely to be an issue (i.e. -9999 values)
CHIRPS<-terra::rast(paste0(CHIRPS_dir,"/",FILES$File[FILES$Year==2010]))
CHIRPS<-terra::as.array(CHIRPS)

CellReferenceAll<-data.table(reshape2::melt(CHIRPS[,,1]))
setnames(CellReferenceAll,c("Var1","Var2"),c("Row","Col"))
CellReferenceAll<-CellReferenceAll[,Index:=1:nrow(CellReferenceAll)
][,is.na:=value==-9999
][,ADMIN:=points$ADMIN
  ][,value:=NULL]


rm(CHIRPS,points,AfricaMap,CHIRPSrast)

# To merge all data across years by site we will need to load all years into memory, one saved transposed year of data is about 4.1gb in memory
# and there are 40 files. So I am splitting each into 10 chunks. Loading 40 x 0.4gb chunks into memory should use about 16gb RAM.

save(CellReferenceAll,file=paste0(SaveDir,"/CellReference.RData"))

Times<-data.table(Date=seq(min(FILES$Date),max(FILES$Date),by="days"))[,N:=1:.N]
save(Times,file=paste0(SaveDir,"/Times.RData"))


for(i in min(FILES$Year):max(FILES$Year)){
  # Display progress
  cat('\r                                                                          ')
  cat('\r',paste0("Loading data: ",i,"/",max(FILES$Year)))
  flush.console()
  
  CHIRPS<-terra::rast(paste0(CHIRPS_dir,"/",FILES$File[FILES$Year==i]))
  CHIRPS<-terra::as.array(CHIRPS)
  
  CHIRPSlist<-asplit(CHIRPS,MARGIN=c(1,2))
  
  # Split list into chunks
  CHIRPSlist<-split(CHIRPSlist, CellReferenceAll$ADMIN)
  
  for(j in 1:length(CHIRPSlist)){
    # Display progress
    cat('\r                                                                                                                                          ')
    cat('\r',paste0("Processing year: ",i,"/",max(FILES$Year), " | ",names(CHIRPSlist)[j]))
    flush.console()
    
    X<-CHIRPSlist[[j]]
    names(X)<-CellReferenceAll[ADMIN==names(CHIRPSlist)[j],Index]
    save(X,file=paste0(SaveDir,"/",names(CHIRPSlist)[j],"-",i))
  }
  
}

# Merge chunks over time
CellReferenceAll<-miceadds::load.Rdata2(file="/CellReference.RData", path=SaveDir)

rm(CHIRPSlist,CHIRPS)
closeAllConnections()
gc()

Countries<-CellReferenceAll[!is.na(ADMIN),unique(ADMIN)]

for(COUNTRY in  Countries){
  
  if(!file.exists(paste0(SaveDir,"/",COUNTRY,".RData"))){
    Files<-list.files(SaveDir,paste0(COUNTRY,"-"))
    
    if(COUNTRY %in% c("Sudan","Republic of the Congo","Guinea")){
      Files<-Files[!grepl("Democratic|South|Bissau|Equatorial",Files)]
    }
    
    X<-lapply(1:length(Files),FUN=function(j){
      # Display progress
      cat('\r                                                                                                                                          ')
      cat(paste0("Loading country: ",COUNTRY," - ",j,"/",length(Files)," - ",Files[j]))
      flush.console()
      
      Y<-load.Rdata2(filename=Files[j],path=SaveDir)
      Y
    })
    
    Y<-pblapply(1:length(X[[1]]),FUN=function(i){
      unlist(lapply(X,"[[",i))
    })
    
    names(Y)<-names(X[[1]])
    
    rm(X)
    gc()
    
    save(Y,file=paste0(SaveDir,"/",COUNTRY,".RData"),compress="gzip",compression_level=6)
    
    rm(Y)
    gc()
  }
  
}

# Delete year chunks
unlink(list.files(SaveDir,"-",full.names = T))

