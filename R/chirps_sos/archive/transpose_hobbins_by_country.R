require(data.table)
require(pbapply)
require(miceadds)
require(terra)
require(sf)
require(sp)

CHIRPS_dir<-"/home/jovyan/common_data/chirps_af/raw"
Data_dir<-"/home/jovyan/common_data/hobbins_ref_et/raw"
Save_dir<-"/home/jovyan/common_data/hobbins_ref_et/intermediate/stacks"
Save_dir2<-"/home/jovyan/common_data/hobbins_ref_et/intermediate/countries"

if(!dir.exists(Save_dir)){
  dir.create(Save_dir,recursive = T)
}

if(!dir.exists(Save_dir2)){
  dir.create(Save_dir2,recursive = T)
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

FILES<-data.table(FileName=list.files(Data_dir,".tif"))
FILES<-FILES[!grepl("_1980",FileName)]

FullPath<-list.files(Data_dir,".tif",full.names = T)
FullPath<-FullPath[!grepl("_1980",FullPath)]

FILES[,Year:=substr(FileName,5,5+3)
      ][,Dekad:=substr(FileName,9,10)
        ][,FullPath:=FullPath]


CHIRPSrast<-terra::rast(list.files(CHIRPS_dir,full.names = T)[1])


for(YEAR in FILES[,unique(Year)]){
  
  SaveFile<-paste0(Save_dir,"/",YEAR,".tif")
  
  if(!file.exists(SaveFile)){
  
    # Display progress
    cat('\r                                                                           ')
    cat('\r',paste0("Processing: ",paste(YEAR,"/",FILES[,max(Year)])))
    flush.console()
    
  
    PET.Year<-terra::rast(FILES[Year==YEAR,FullPath])
    PET.Year<-terra::resample(PET.Year,CHIRPSrast)
    
    terra::writeRaster(PET.Year,file=SaveFile)
  }

}

FILES<-list.files(Save_dir,full.names = T)

PET.Year<-terra::as.array(terra::rast(FILES[10]))

Geom<-data.table(terra::crds(CHIRPSrast))

colnames(PET.Year)<-Geom[,unique(x)]
rownames(PET.Year)<-Geom[,unique(y)]

# Load map of Africa
AfricaMap<-rworldmap::getMap(resolution = "high")
AfricaMap<-AfricaMap[AfricaMap$REGION=="Africa"&!is.na(AfricaMap$REGION),]
AfricaMap<-sf::st_as_sf(AfricaMap)
AfricaMap<-terra::vect(AfricaMap)
AfricaMap<-terra::project(AfricaMap,CHIRPSrast)
AfricaMap<-terra::rasterize(AfricaMap,CHIRPSrast,field="ADMIN")

# Convert to points
points<-terra::as.points(AfricaMap,values=T,na.rm=F)

CellReferenceAll<-data.table(reshape2::melt(PET.Year[,,1]))
CellReferenceAll<-CellReferenceAll[,Index:=1:nrow(CellReferenceAll)]
setnames(CellReferenceAll,c("Var1","Var2"),c("Row","Col"))

CellReferenceAll<-CellReferenceAll[,is.na:=value==-9999
                                   ][,ADMIN:=points$ADMIN
                                     ][,value:=NULL]

XY<-data.table(terra::geom(terra::as.points(CHIRPSrast, values=TRUE, na.rm=F)))[,list(x,y)]
CellReferenceAll<-cbind(CellReferenceAll,XY)


save(CellReferenceAll,file=paste0(Save_dir2,"/CellReference.RData"))
Times<-data.table(Year=rep(1981:2021,each=36),Dekad=rep(1:36,length(1981:2021)))[,N:=1:.N]
save(Times,file=paste0(Save_dir2,"/Times.RData"))

Files<-list.files(Save_dir)

# Turn raster into 
for(i in 1:length(Files)){
  
    cat('\r                                                                                                                                          ')
    cat('\r',paste0("Loading data: ",i,"/",length(Files))," - ",Files[i])
    flush.console()
    
  
    PET<-terra::as.array(terra::rast(paste0(Save_dir,"/",Files[i])))
    colnames(PET)<-Geom[,unique(x)]
    rownames(PET)<-Geom[,unique(y)]
    
    # Split array into a list, each entry in the list is a non-NA raster cell
    PETlist<-asplit(PET,MARGIN=c(1,2))

    # Split list into chunks
    PETlist<-split(PETlist, CellReferenceAll$ADMIN)
    
    for(j in 1:length(PETlist)){
      
      # Display progress
      cat('\r                                                                                                                                          ')
      cat('\r',paste0("Processing file: ",i,"/",length(Files))," - ",Files[i], " | ",names(PETlist)[j])
      flush.console()
      
      X<-PETlist[[j]]
      names(X)<-CellReferenceAll[ADMIN==names(PETlist)[j],Index]
      save(X,file=paste0(Save_dir2,"/",names(PETlist)[j],"-",Files[i]))
    }
  }


# Merge chunks over time
CellReferenceAll<-miceadds::load.Rdata2(file="/CellReference.RData", path=Save_dir2)


rm(PETlist)
closeAllConnections()
gc()

Countries<-CellReferenceAll[!is.na(ADMIN),unique(ADMIN)]

for(COUNTRY in  Countries){
  
  if(!file.exists(paste0(Save_dir2,"/",COUNTRY,".RData"))){
    Files<-list.files(Save_dir2,paste0(COUNTRY,"-"))
    
    if(COUNTRY %in% c("Sudan","Republic of the Congo","Guinea")){
      Files<-Files[!grepl("Democratic|South|Bissau|Equatorial",Files)]
    }
    
    X<-lapply(1:length(Files),FUN=function(j){
      # Display progress
      cat('\r                                                                                                          ')
      cat(paste0("Loading country: ",COUNTRY," - ",j,"/",length(Files)," - ",Files[j]))
      flush.console()
      
      
      Y<-load.Rdata2(filename=Files[j],path=Save_dir2)
      Y
    })
    
    Y<-pblapply(1:length(X[[1]]),FUN=function(i){
      unlist(lapply(X,"[[",i))
    })
    
    names(Y)<-names(X[[1]])
    
    save(Y,file=paste0(Save_dir2,"/",COUNTRY,".RData"),compress="gzip",compression_level=6)
    
    rm(X)
    rm(Y)
    gc()
  }
  
  
}


# Delete year chunks
unlink(list.files(Save_dir2,"-",full.names = T))



