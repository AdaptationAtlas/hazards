require(ggplot2)
require(data.table)
require(miceadds)
require(sp)
require(pbapply)
require(raster)
require(rworldmap)
require(rworldxtra)

# Version
version<-1


DataDir<-"C:/Datasets"

# Set directories
CHIRPSraw_Dir<-paste0(DataDir,"/chirps_af/raw")
CHIRPS_Dekad_Dir<-paste0(DataDir,"/chirps_af/intermediate/array_countries_dekad")
SOS_Dir<-paste0(DataDir,"/atlas_SOS/intermediate/v",version)
list.dirs(SOS_Dir)
SOS_Dir<-list.dirs(SOS_Dir)[2]

# Load map of Africa
CHIRPSrast<-terra::rast(list.files(CHIRPSraw_Dir,full.names = T)[1])

AfricaMap<-rworldmap::getMap(resolution = "high")
AfricaMap<-AfricaMap[AfricaMap$REGION=="Africa"&!is.na(AfricaMap$REGION),]
AfricaMap<-sf::st_as_sf(AfricaMap)
AfricaMap<-terra::vect(AfricaMap)
AfricaMap<-terra::project(AfricaMap,CHIRPSrast)

Countries<-as.character(AfricaMap$ADMIN)
NotInHobbins<-c("Cape Verde","Mauritius","Saint Helena","Seychelles")
Issue<-c("Djibouti","Western Sahara")
Countries<-Countries[!Countries %in% c(NotInHobbins,Issue)]

# Wrangle SOS chunks ####
SOSFiles<-data.table(File=list.files(SOS_Dir,".RData",full.names = T),
                     Country=gsub("[.]RData","",tstrsplit(list.files(SaveDir,".RData",full.names = F),"_",keep=2)[[1]]))

AnalysisDir<-paste0(SaveDir,"/explore_results")
if(!dir.exists(AnalysisDir)){
  dir.create(AnalysisDir)
}

File<-paste0(AnalysisDir,"/SOSlist.RData")

# Join LTAvg data across countries
if(!file.exists(File)){
  SOSlist<-pblapply(1:nrow(SOSFiles),FUN=function(i){
    load(file=SOSFiles[i,File])
    list(LTAvg2=SOS_Data$LTAvg_SOS2,
         LTAvg3=SOS_Data$LTAvg_SOS3,
         #Seasonal2=SOS_Data$Seasonal_SOS2,
         #Seasonal3=SOS_Data$Seasonal_SOS3,
         LTAVg_Data=SOS_Data$LTAVg_Data,
         LTAVg_Summary=SOS_Data$LTAVg_Summary)
  })
  
  names(SOSlist)<-SOSFiles$Country
  
  save(SOSlist,file=File)
}else{
  load(File)
}


gc()


#***************************************#
# Calculate only from LTAVg_Summary ####
MinLength<-1 # Minimum length of season allowed
Min.Rain<-0 # Minimum seasonal rainfall allowed
SeasonGapT<-3 # if the gap between seasons is <=SeasonGapT then it is considered as a split season (i.e. bimodal)

PlotDir1<-paste0(SaveDir,"/Plots/ML",MinLength,"-MR",Min.Rain,"-SG",SeasonGapT)

# The code will modify EOS, SOS and LGP according to SeasonGapT and attempt to tidy seasonal organization between S1 and S2

Statistic<-"mode" 
Dataset<-"Seasonal"

if(Dataset=="LT"){
  PlotDir<-paste0(PlotDir1,"/LT")
}else{
  PlotDir<-paste0(PlotDir1,"/Seasonal-",Statistic)
}

if(!dir.exists(PlotDir)){
  dir.create(PlotDir,recursive=T)
}

  
for(COUNTRY in Countries){
  print(COUNTRY)
  
  Map.Subset<-AfricaMap[AfricaMap$ADMIN==COUNTRY,]
  
  if(Dataset=="LT"){
   SOSData1<-data.table::copy(SOSlist[[COUNTRY]]$LTAVg_Summary)
  }else{
      SOSData1<-data.table::copy(SOSlist[[COUNTRY]]$LTAvg2)
    setnames(SOSData1,c("Dekad.Season","Tot.Rain.mean"),c("Seq","Tot.Rain"))
    SOSData1[,Seasons:=length(unique(Seq)),by=Index]
  
    if(Statistic=="mean"){
      SOSData1[,SOS:=round(SOS.mean,0)][SOS==0,SOS:=36]
      SOSData1[,EOS:=round(EOS.mean,0)][EOS==0,EOS:=36]
      SOSData1[,LGP:=round(LGP.mean,0)]
    }
    
    if(Statistic=="median"){
      SOSData1[,SOS:=round(SOS.median,0)][SOS==0,SOS:=36]
      SOSData1[,EOS:=round(EOS.median,0)][EOS==0,EOS:=36]
      SOSData1[,LGP:=round(LGP.median,0)]
    }
    
    if(Statistic=="mode"){
      SOSData1[,SOS:=round(SOS.mode,0)][SOS==0,SOS:=36]
      SOSData1[,EOS:=round(EOS.mode,0)][EOS==0,EOS:=36]
      SOSData1[,LGP:=round(LGP.mode,0)]
    }
  }

  # NOTE - here seasons with insufficient length or rainfall are removed. This could be done after merging of bimodal seasons?
  SOSData1<-SOSData1[LGP>=MinLength & Tot.Rain>Min.Rain]

  if(Dataset=="LT"){
   suppressWarnings(SOSData1[,Seq:=stringi::stri_replace_all_regex(Seq,
                                                                  pattern=unique(Seq[!is.na(Seq)]),
                                                                   replacement=1:length(unique(Seq[!is.na(Seq)]))),by=Index])
  }
  
  # Load country raster
  CHIRPS_Ref<-terra::rast(paste0(CHIRPS_Dekad_Dir,"/",COUNTRY,".tif"))[[1]]
  CHIRPS_Ref<-data.table(terra::as.data.frame(CHIRPS_Ref,xy=T,cells=T))[,1:3]
  setnames(CHIRPS_Ref,"cell","Index")
  
  SOSData1<-merge(CHIRPS_Ref,SOSData1,by="Index",all.x=F)
  
  SOSData1<-dcast(SOSData1,Index+x+y+Seasons~Seq,value.var = c("SOS","EOS","LGP","Tot.Rain"))
  if(!"SOS_2" %in% colnames(SOSData1)){
    SOSData1[,c("SOS_2","EOS_2","LGP_2","Tot.Rain_2","Tot.ETo_2","Balance_2"):=NA]
  }
  
  # Calculate distance between seasons
  SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist12:=round(CicularDist(EOS_1*10,SOS_2*10)/10,0),by=Index]
  SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist21:=round(CicularDist(EOS_2*10,SOS_1*10)/10,0),by=Index]
  
  # Determine the order of seasons by the length of the distance between them, this is used to determine how bimodal
  # seasons are merged
  SOSData1[,MinDist:=paste(which(c(Dist12,Dist21)<=SeasonGapT),collapse="-"),by=Index][MinDist=="",MinDist:=NA]
  SOSData1[,MinDistN:=sum(c(Dist12,Dist21)<=SeasonGapT),by=Index]
  
  SOSData<-data.table::copy(SOSData1)
  
  # Note: Consider using S1, S2, etc. to record modified values
  SOSData[,S1:=SOS_1][,S2:=SOS_2][,E1:=EOS_1][,E2:=EOS_2][,LGP1:=LGP_1][,LGP2:=LGP_2]
  
  SOSData[MinDistN==2,E1:=if(S1==1){36}else{S1-1},by=Index]
  
  SOSData[MinDistN==1 & MinDist==1,E1:=E2]
  
  SOSData[MinDistN==1 & MinDist==2,S1:=S2]
  
  SOSData[MinDistN %in% c(1,2),LGP1:=LGP1+LGP2
  ][MinDistN %in% c(1,2),LGP2:=NA
  ][MinDistN %in% c(1,2),E2:=NA
  ][MinDistN %in% c(1,2),S2:=NA]
  
  SOSData[,Seasons:=sum(c(!is.na(S1),!is.na(S2)),na.rm=T),by=Index]
  
  if(SOSData[,!any(Seasons==2)]){
    SEASONS<-list(S1=OrderDekadSeq(unique(round(SOSData$S1))),S2=NULL)
    
    SOSData[,S1:=round(S1,0)][,S2:=NA][S1==0,S1:=36]
    
    SOSData[,Swettest:=S1][!is.na(S1),Seasons:=1]
  }else{
    
    if(SOSData[,!any(Seasons==1)]){
      SEASONS<-list(S1=OrderDekadSeq(unique(round(SOSData$S2))),S2=NULL)
      
      SOSData[,S1:=round(S1,0)][,S2:=NA][S1==0,S1:=36]
      
      SOSData[,Swettest:=S1][!is.na(S1),Seasons:=1]
    }else{
      
      
      # Use multimode function to determine seasons from the distribution of SOS values (from all cells and both seasons combined)
      # This approach doesn't work well for some countries and seasons are specified manually
      if(!COUNTRY %in% c("Kenya","South Africa","Mozambique")){
        Modes<-suppressWarnings(multimode::locmodes(data=c(SOSData$S1,SOSData$S2),mod0=2,display=TRUE))
        Modes<-round(Modes$locations,0)
        S1<-Modes[1]:Modes[2]
        S2<-(Modes[2]+1):Modes[3]
        Left<-if(Modes[1]==1){NULL}else{1:(Modes[1]-1)}
        Right<-if(Modes[3]==36){NULL}else{(Modes[3]+1):36}
        Rem<-c(Right,Left)
        
        if(length(Rem) %% 2 == 0){
          S2<-c(S2,Rem[1:(length(Rem)/2)])
          S1<-c(Rem[(1+(length(Rem)/2)):length(Rem)],S1)
        }else{
          S2<-c(S2,Rem[1:(floor(length(Rem)/2))])
          S1<-c(Rem[ceiling(((length(Rem)/2))):length(Rem)],S1)
        }
        
        SEASONS<-list(S1=S1,S2=S2)
      }
      
      if(COUNTRY=="Kenya"){
        SEASONS<-list(S1=c(35:36,1:15),S2=16:34)
      }
      if(COUNTRY=="South Africa"){
        SEASONS<-list(S1=8:20,S2=c(21:36,1:7))
      }
      if(COUNTRY=="Mozambique"){
        SEASONS<-list(S1=1:36,S2=NULL)
      }
      
      SOSData[,X1:=as.numeric(NA)
      ][,X2:=as.numeric(NA)
      ][S1 %in% SEASONS$S1,X1:=S1
      ][S1 %in% SEASONS$S2,X2:=S1
      ][S2 %in% SEASONS$S1,X1:=S2
      ][S2 %in% SEASONS$S2,X2:=S2
      ][S1 %in% SEASONS$S1 & S2 %in% SEASONS$S1,X1:=S1
      ][S1 %in% SEASONS$S1 & S2 %in% SEASONS$S1,X2:=S2 
      ][S1 %in% SEASONS$S2 & S2 %in% SEASONS$S2,X1:=S1
      ][S1 %in% SEASONS$S2 & S2 %in% SEASONS$S2,X2:=S2]
      
      # Assign LGP to reorganised seasons
      SOSData[X1==S1,Y1:=LGP1
      ][X1==S2,Y1:=LGP2
      ][X2==S1,Y2:=LGP1
      ][X2==S2,Y2:=LGP2]
      
      SOSData[,S1:=X1][,S2:=X2][,X1:=NULL][,X2:=NULL][,LGP1:=Y1][,LGP2:=Y2][,Y1:=NULL][,Y2:=NULL]
      
      SOSData[,Swettest:=SOS_1
      ][!is.na(SOS_2) & Tot.Rain_1>Tot.Rain_2,Swettest:=SOS_2
      ][is.na(Swettest),Swettest:=SOS_2]
      
    }}
  
  SOSData[,Seasons:=sum(c(!is.na(S1),!is.na(S2)))*1,by=Index]
  
  # Manually tidy up countries which are largely unimodal but have a few pixels that split the seasons.
  # Note: Check that correct LGP field is being used here
  if(COUNTRY %in% c("Nigeria","Burundi","Madagascar","Mozambique")){
    SOSData[is.na(S1) & !is.na(S2),EOS_1:=EOS_2
    ][is.na(S1) & !is.na(S2),EOS_2:=NA
    ][is.na(S1) & !is.na(S2),LGP1 :=LGP2
    ][is.na(S1) & !is.na(S2),LGP2:=NA
    ][is.na(S1) & !is.na(S2),Tot.Rain_1:=Tot.Rain_2
    ][is.na(S1) & !is.na(S2),Tot.Rain_2:=NA
    ][is.na(S1) & !is.na(S2),Tot.ETo_1:=Tot.ETo_2
    ][is.na(S1) & !is.na(S2),Tot.ETo_2:=NA
    ][is.na(S1) & !is.na(S2),Balance_1:=Balance_2
    ][is.na(S1) & !is.na(S2),Balance_2:=NA
    ][is.na(S1) & !is.na(S2),S1:=S2
    ][S1==S2,S2:=NA]
  }
  
  # Classify systems
  SOSData[,LGPsum:=sum(c(LGP1,LGP2),na.rm=T),by=Index
  ][Seasons==2,System:="2 WS"
  ][Seasons==1 & MinDistN %in% c(1,2),System:="1 WS-Split"
  ][Seasons==1 & !MinDistN %in% c(1,2),System:="1 WS"
  ][LGPsum>=30,System:=paste0(System,"-LGP>=30")
  ][LGPsum<=6,System:=paste0(System,"-LGP<=6")]
  
  if(!dir.exists(paste0(PlotDir,"/Data"))){
    dir.create(paste0(PlotDir,"/Data"))
  }
  Colnames<-colnames(SOSData)[!grepl("_3|_4|_5|_6",colnames(SOSData))]
  SOSData<-SOSData[,..Colnames]
  save(SOSData,file=paste0(PlotDir,"/Data/",COUNTRY,".RData"))
  
  # prepare coordinates, data, and proj4string
  X <- terra::vect(SOSData,geom=c("x", "y"),crs=terra::crs(CHIRPSrast))
  
  BaseRaster<-terra::rast(ncol=SOSData[,length(unique(x))], 
                          nrow=SOSData[,length(unique(y))], 
                          xmin=SOSData[,min(x)], 
                          xmax=SOSData[,max(x)], 
                          ymin=SOSData[,min(y)],
                          ymax=SOSData[,max(y)])
  
  # SOS Maps ####
  Y1<-round(terra::rasterize(X,BaseRaster,field="S1"),0)
  if(SOSData[!is.na(S2),.N==0]){
    Y2<-Y1
    Y2[]<-NA
  }else{
    Y2<-round(terra::rasterize(X,BaseRaster,field="S2"),0)
  }
  Ymax<-round(terra::rasterize(X,BaseRaster,field="Swettest"),0)
  
  LGP_1<-round(terra::rasterize(X,BaseRaster,field="LGP1"),0)
  
  if(!2 %in% SOSData$Seasons){
    LGP_2<-LGP_1
    LGP_2[]<-NA
  }else{
    LGP_2<-round(terra::rasterize(X,BaseRaster,field="LGP2"),0)
  }
  
  
  SOS<-c(Y1,Y2,Ymax)
  LGP<-c(LGP_1,LGP_2,sum(c(LGP_1,LGP_2),na.rm = T))
  
  names(SOS)<-c("SOS_1",
                "SOS_2",
                "SOS_Wettest")
  
  names(LGP)<-c("LGP_1",
                "LGP_2",
                "LGP_sum")
  
  SOS2<-raster::stack(SOS)
  LGP2<-raster::stack(LGP)
  
  coords <- xyFromCell(SOS2, seq_len(ncell(SOS2)))
  Data <- raster::stack(as.data.frame(raster::getValues(SOS2)))
  names(SOS2) <- names(SOS)
  
  Data <- cbind(coords, Data)
  Data$values<-factor(Data$values,levels=1:36)
  
  SOSmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=3) +
    scale_fill_manual(
      values=pals::kovesi.cyclic_mygbm_30_95_c78_s25(n=36),
      na.value = "transparent",
      drop = FALSE,
      na.translate = F
    )+
    theme_bw()+
    labs(fill="SOS")+
    theme(legend.position = "bottom",
          axis.title = element_blank())+
    guides(fill=guide_legend(nrow=3))+
    coord_equal()
  
  coords <- xyFromCell(LGP2, seq_len(ncell(LGP2)))
  Data <- raster::stack(as.data.frame(raster::getValues(LGP2)))
  names(LGP2) <- names(LGP)
  
  Data <- cbind(coords, Data)

  LGPmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=3) +
    viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
    theme_bw()+
    labs(fill="LGP (dekads)")+
    theme(legend.position = "bottom",
          axis.title = element_blank())+
    coord_equal()
  
  # Seasonality Map ####
  
  Ysystem<-terra::rasterize(X,BaseRaster,field="System")
  names(Ysystem)<-"System"
  coords <- xyFromCell(Ysystem, seq_len(ncell(Ysystem)))
  Data <-values(Ysystem)
  names(Ysystem) <- names(Ysystem)
  
  Levels<-data.table(value=sort(unique(Ysystem[!is.na(Ysystem)])),description=levels(Ysystem)[[1]])
  
  Data <- data.frame(cbind(coords, Data))
  Data$System<-Levels$description.category[match(Data$System,Levels$value)]
  
  Ysystem<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = System)) +
    viridis::scale_fill_viridis(option="turbo",discrete = T,na.value="transparent",na.translate = F)+
    theme_bw()+
    labs(fill="System")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.direction = "vertical")+
    guides(fill=guide_legend(ncol=1))+
    coord_equal()
  
  # No Seasons Map ####
  Yseasons<-terra::rasterize(X,BaseRaster,field="Seasons")
  names(Yseasons)<-"Seasons"
  coords <- xyFromCell(Yseasons, seq_len(ncell(Yseasons)))
  Data <-values(Yseasons)
  names(Yseasons) <- names(Yseasons)
  
  Data <- data.frame(cbind(coords, Seasons=Data))
  
  Yseasons<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = as.factor(Seasons))) +
    viridis::scale_fill_viridis(option="viridis",discrete = T,na.value="transparent",na.translate = F,direction=-1)+
    theme_bw()+
    labs(fill="Seasons")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.direction = "vertical")+
    guides(fill=guide_legend(ncol=1))+
    coord_equal()
  
  
  g1<-gridExtra::arrangeGrob(grobs=list(SOSmap,LGPmap),layout_matrix=matrix(rep(1:2,3),nrow = 2))
  g2<-gridExtra::arrangeGrob(grobs=list(Yseasons,g1,Ysystem),layout_matrix=matrix(c(1,rep(2,3),3),nrow = 1))
  
  ggsave(filename = paste0(COUNTRY,".png"),
         plot = g2,
         path = PlotDir,
         width= 200,
         height = 100,
         units = "mm",
         scale = 1.8,
         device=png,
         dpi = 600,
         bg="white")
  
}

# Plot all countries ####
Files<-list.files(paste0(PlotDir,"/Data"),".RData")

Data<-miceadds::load.Rdata2(filename=Files[1],path=paste0(PlotDir,"/Data"))

SOSData<-data.table::rbindlist(lapply(Files,FUN=function(X){miceadds::load.Rdata2(filename=X,path=paste0(PlotDir,"/Data"))}),use.names = T)
SOSData[,Seasons:=sum(!is.na(SOS_1),!is.na(SOS_2))]

# prepare coordinates, data, and proj4string
coords <- data.frame(SOSData[ , c("x", "y")])   # coordinates
crs    <- CRS("+init=epsg:4326") # proj4string of coords

# make the SpatialPointsDataFrame object
X <- terra::vect(SOSData,geom=c("x", "y"),crs=terra::crs(CHIRPSrast))

BaseRaster<-terra::rast(ncol=SOSData[,length(unique(x))], 
                        nrow=SOSData[,length(unique(y))], 
                        xmin=SOSData[,min(x)], 
                        xmax=SOSData[,max(x)], 
                        ymin=SOSData[,min(y)],
                        ymax=SOSData[,max(y)])

# SOS Maps ####
# Note here we are using the unmodified values - no effort has been made to seasons into season 1/2
# Season can be merged according to rules
Y1<-round(terra::rasterize(X,BaseRaster,field="SOS_1"),0)
names(Y1)<-"SOS_1"
terra::writeRaster(Y1,paste0(PlotDir,"/Data/Africa_SOS1.tif"),overwrite=T)

Y2<-round(terra::rasterize(X,BaseRaster,field="SOS_2"),0)
names(Y2)<-"SOS_2"
terra::writeRaster(Y2,paste0(PlotDir,"/Data/Africa_SOS2.tif"),overwrite=T)

Ymax<-round(terra::rasterize(X,BaseRaster,field="Swettest"),0)
names(Ymax)<-"SOS_Wettest"
terra::writeRaster(Ymax,paste0(PlotDir,"/Data/Africa_SOSmaxrain.tif"),overwrite=T)

Yseasons<-terra::rasterize(X,BaseRaster,field="Seasons")
names(Yseasons)<-"Seasons"
terra::writeRaster(Yseasons,paste0(PlotDir,"/Data/Africa_Seasons.tif"),overwrite=T)

# LGP Maps
LGP1<-round(terra::rasterize(X,BaseRaster,field="LGP_1"),0)
names(LGP1)<-"LGP_1"
terra::writeRaster(LGP1,paste0(PlotDir,"/Data/Africa_LGP1.tif"),overwrite=T)

LGP2<-round(terra::rasterize(X,BaseRaster,field="LGP_2"),0)
names(LGP2)<-"LGP_2"
terra::writeRaster(LGP2,paste0(PlotDir,"/Data/Africa_LGP2.tif"),overwrite=T)

# EOS Maps
EOS1<-round(terra::rasterize(X,BaseRaster,field="EOS_1"),0)
names(EOS1)<-"EOS_1"
terra::writeRaster(EOS1,paste0(PlotDir,"/Data/Africa_EOS1.tif"),overwrite=T)

EOS2<-round(terra::rasterize(X,BaseRaster,field="EOS_2"),0)
names(EOS2)<-"EOS_2"
terra::writeRaster(EOS2,paste0(PlotDir,"/Data/Africa_EOS2.tif"),overwrite=T)

# Stack MAPS

SOS<-c(Y1,Y2,Ymax)
LGP<-c(LGP1,LGP2,sum(c(LGP1,LGP2),na.rm = T))

names(SOS)<-c("SOS_1",
              "SOS_2",
              "SOS_Wettest")

names(LGP)<-c("LGP_1",
              "LGP_2",
              "LGP_sum")

SOS2<-raster::stack(SOS)
LGP2<-raster::stack(LGP)

coords <- xyFromCell(SOS2, seq_len(ncell(SOS2)))
Data <- raster::stack(as.data.frame(raster::getValues(SOS2)))
names(SOS2) <- names(SOS)

Data <- cbind(coords, Data)
Data$values<-factor(Data$values,levels=1:36)

SOSmap<-ggplot(Data) + 
  geom_tile(aes(x, y, fill = values)) +
  facet_wrap(~ ind,ncol=3) +
  scale_fill_manual(
    values=pals::kovesi.cyclic_mygbm_30_95_c78_s25(n=36),
    na.value = "transparent",
    drop = FALSE,
    na.translate = F
  )+
  theme_bw()+
  labs(fill="SOS")+
  theme(legend.position = "bottom",
        axis.title = element_blank())+
  guides(fill=guide_legend(nrow=3))+
  coord_equal()

ggsave(filename = "Africa - SOS.png",
       plot = SOSmap,
       path = PlotDir,
       width= 200,
       height = 100,
       units = "mm",
       scale = 1.2,
       dpi = 600,
       device = png,
       bg="white")

# LGP ####

coords <- xyFromCell(LGP2, seq_len(ncell(LGP2)))
Data <- raster::stack(as.data.frame(raster::getValues(LGP2)))
names(LGP2) <- names(LGP)

Data <- cbind(coords, Data)

LGPmap<-ggplot(Data) + 
  geom_tile(aes(x, y, fill = values)) +
  facet_wrap(~ ind,ncol=3) +
  viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
  theme_bw()+
  labs(fill="LGP (dekads)")+
  theme(legend.position = "bottom",
        axis.title = element_blank())+
  coord_equal()

ggsave(filename = "Africa - LGP.png",
       plot = LGPmap,
       path = PlotDir,
       width= 200,
       height = 100,
       units = "mm",
       scale = 1.2,
       dpi = 600,
       device = png,
       bg="white")

# Seasonality Map ####

Ysystem<-terra::rasterize(X,BaseRaster,field="System")
names(Ysystem)<-"System"
coords <- xyFromCell(Ysystem, seq_len(ncell(Ysystem)))
Data <-values(Ysystem)
names(Ysystem) <- names(Ysystem)

Levels<-data.table(value=sort(unique(Ysystem[!is.na(Ysystem)])),description=levels(Ysystem)[[1]])

Data <- data.frame(cbind(coords, Data))
Data$System<-Levels$description.category[match(Data$System,Levels$value)]

Ysystem<-ggplot(Data) + 
  geom_tile(aes(x, y, fill = System)) +
  viridis::scale_fill_viridis(option="turbo",discrete = T,na.value="transparent",na.translate = F)+
  theme_bw()+
  labs(fill="System")+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        legend.direction = "vertical")+
  guides(fill=guide_legend(ncol=1))+
  coord_equal()

Ysystem2<-Ysystem+theme(legend.position = "bottom",
              axis.title = element_blank(),
              legend.direction = "horizontal")+
  guides(fill=guide_legend(nrow=2))

ggsave(filename = "Africa - systems.png",
       plot = Ysystem2,
       path = PlotDir,
       width= 120,
       height = 120,
       units = "mm",
       scale = 1.6,
       dpi = 600,
       device = png,
       bg="white")

g1<-gridExtra::arrangeGrob(grobs=list(SOSmap,LGPmap),layout_matrix=matrix(rep(1:2,3),nrow = 2))
g2<-gridExtra::arrangeGrob(grobs=list(g1,Ysystem),layout_matrix=matrix(c(rep(1,3),2),nrow = 1))


ggsave(filename = "Africa - Combined.png",
       plot = g2,
       path = PlotDir,
       width= 200,
       height = 100,
       units = "mm",
       scale = 1.8,
       dpi = 600,
       device = png,
       bg="white")







#*****************************************************#
#####
#####
#####
# Calculate across the entire time series of data ####
for(COUNTRY in Countries){
  
  print(COUNTRY)  
  
  Map.Subset<-AfricaMap[grep(COUNTRY,AfricaMap$ADMIN),]
  
  SOSData<-data.table::copy(SOSlist[[grep(COUNTRY,SOSFiles)]]$LTAvg2)
  SOSData.S1<-merge(CHIRPS_Ref,SOSData[Dekad.Season==1],by="Index",all.x=T)
  SOSData.S2<-merge(CHIRPS_Ref,SOSData[Dekad.Season==2],by="Index",all.x=T)
  
  SOSData<-data.table::copy(SOSlist[[grep(COUNTRY,SOSFiles)]]$LTAvg2)
  SOSData<-dcast(SOSData,Index~Dekad.Season,value.var = c("SOS.mean","SOS.sd","Tot.Rain.mean","LGP.mean"))
  
  if(is.null(SOSData$SOS.mean_2)){
    Seasons<-list(S1=OrderDekadSeq(unique(round(SOSData$SOS.mean_1))),S2=NULL)
    
    SOSData[,S1:=round(SOS.mean_1,0)][,S2:=NA][S1==0,S1:=36]
    
    SOSData[,Swettest:=S1][!is.na(S1),Seasons:=1]
    
    # Assign LGP to reorganised seasons
    SOSData[,LGP1:=LGP.mean_1][,LGP2:=NA]
    
  }else{
    # Seasons<- FindSeasonFun(Data=c(SOSData$SOS.mean_1,SOSData$SOS.mean_2))
    
    Modes<-multimode::locmodes(data=c(SOSData$SOS.mean_1,SOSData$SOS.mean_2),mod0=2,display=TRUE)
    Modes<-round(Modes$locations,0)
    S1<-Modes[1]:Modes[2]
    S2<-(Modes[2]+1):Modes[3]
    Left<-if(Modes[1]==1){NULL}else{1:(Modes[1]-1)}
    Right<-if(Modes[3]==36){NULL}else{(Modes[3]+1):36}
    Rem<-c(Right,Left)
    
    if(length(Rem) %% 2){
      S2<-c(S2,Rem[1:(length(Rem)/2)])
      S1<-c(Rem[(1+(length(Rem)/2)):length(Rem)],S1)
    }else{
      S2<-c(S2,Rem[1:(floor(length(Rem)/2))])
      S1<-c(Rem[ceiling(((length(Rem)/2))):length(Rem)],S1)
    }
    
    Seasons<-list(S1=S1,S2=S2)
    
    SOSData[,SOS.mean_1:=round(SOS.mean_1,0)
    ][,SOS.mean_2:=round(SOS.mean_2,0)
    ][SOS.mean_1==0,SOS.mean_1:=36
    ][SOS.mean_2==0,SOS.mean_2:=36]
    
    SOSData[SOS.mean_1 %in% Seasons$S1,S1:=SOS.mean_1
    ][SOS.mean_1 %in% Seasons$S2,S2:=SOS.mean_1
    ][SOS.mean_2 %in% Seasons$S1,S1:=SOS.mean_2
    ][SOS.mean_2 %in% Seasons$S2,S2:=SOS.mean_2]
    
    SOSData[SOS.mean_1 %in% Seasons$S2 & SOS.mean_2 %in% Seasons$S2,S2:=min(S1,S2),by=Index]
    SOSData[SOS.mean_1 %in% Seasons$S2 & SOS.mean_2 %in% Seasons$S2,S1:=NA,by=Index]
    
    SOSData[SOS.mean_1 %in% Seasons$S1 & SOS.mean_2 %in% Seasons$S1,S1:=min(S1,S2),by=Index]
    SOSData[SOS.mean_1 %in% Seasons$S1 & SOS.mean_2 %in% Seasons$S1,S2:=NA,by=Index]
    
    SOSData[,Swettest:=S1
    ][!is.na(SOS.mean_2) & Tot.Rain.mean_2>Tot.Rain.mean_1,Swettest:=S2
    ][is.na(Swettest),Swettest:=S2
    ]
    
    # Assign LGP to reorganised seasons
    SOSData[S1==SOS.mean_1,LGP1:=LGP.mean_1
    ][S1==SOS.mean_2,LGP1:=LGP.mean_2
    ][S2==SOS.mean_1,LGP2:=LGP.mean_1
    ][S2==SOS.mean_2,LGP2:=LGP.mean_2]
  }
  
  
  SOSData[,Seasons:=sum(c(!is.na(S1),!is.na(S2)))*10,by=Index]
  
  SOSData.S1<-merge(CHIRPS_Ref,SOSData,by="Index",all.x=F)
  
  X <- terra::vect(SOSData.S1,geom=c("Col", "Row"),crs=terra::crs(CHIRPSrast))
  
  
  BaseRaster<-terra::rast(ncol=SOSData.S1[,length(unique(Col))], 
                          nrow=SOSData.S1[,length(unique(Row))], 
                          xmin=SOSData.S1[,min(Col)], 
                          xmax=SOSData.S1[,max(Col)], 
                          ymin=SOSData.S1[,min(Row)],
                          ymax=SOSData.S1[,max(Row)])
  
  # SOS Maps ####
  Y1<-round(terra::rasterize(X,BaseRaster,field="S1"),0)
  if(is.null(SOSData$SOS.mean_2)){
    Y2<-Y1
    Y2[]<-NA
  }else{
    Y2<-round(terra::rasterize(X,BaseRaster,field="S2"),0)
  }
  Ymax<-round(terra::rasterize(X,BaseRaster,field="Swettest"),0)
  Yseasons<-round(terra::rasterize(X,BaseRaster,field="Seasons"),0)
  
  Ystack<-terra::crop(c(Y1,Y2,Ymax,Yseasons),Map.Subset)
  
  names(Ystack)<-c(paste0(COUNTRY,"_SOS_1"),
                   paste0(COUNTRY,"_SOS_2"),
                   paste0(COUNTRY,"_SOS_Wettest"),
                   paste0(COUNTRY,"_Seasons"))
  
  Ystack2<-raster::stack(Ystack)
  
  coords <- xyFromCell(Ystack2, seq_len(ncell(Ystack2)))
  Data <- raster::stack(as.data.frame(raster::getValues(Ystack2)))
  names(Ystack2) <- names(Ystack)
  
  Data <- cbind(coords, Data)
  Data$values<-factor(Data$values,levels=1:36)
  
  SOSmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=4) +
    #viridis::scale_fill_viridis(option="turbo",discrete = T,drop=F) +
    scale_fill_manual(
      values=pals::kovesi.cyclic_mygbm_30_95_c78_s25(n=36),
      na.value = "transparent",
      drop = FALSE,
      na.translate = F
    )+
    coord_equal()+
    theme_bw()+
    labs(fill="SOS")+
    theme(legend.position = "bottom",
          axis.title = element_blank())+
    guides(fill=guide_legend(nrow=3))
  
  # LGP Maps ####
  LGP1<-round(terra::rasterize(X,BaseRaster,field="LGP1"),0)
  
  if(is.null(SOSData$SOS.mean_2)){
    LGP2<-LGP1
    LGP2[]<-NA
  }else{
    LGP2<-round(terra::rasterize(X,BaseRaster,field="LGP2"),0)
  }
  
  
  LGPstack<-terra::crop(c(LGP1,LGP2),Map.Subset)
  
  names(LGPstack)<-c(paste0(COUNTRY,"_LGP_1"),
                     paste0(COUNTRY,"_LGP_2"))
  
  LGPstack<-raster::stack(LGPstack)
  
  coords <- xyFromCell(LGPstack, seq_len(ncell(LGPstack)))
  Data <- raster::stack(as.data.frame(raster::getValues(LGPstack)))
  names(LGPstack) <- names(LGPstack)
  
  Data <- cbind(coords, Data)
  
  LGPmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=2) +
    viridis::scale_fill_viridis(option="turbo",discrete = F,na.value="transparent")+
    coord_equal()+
    theme_bw()+
    labs(fill="LGP")+
    theme(legend.position = "bottom",
          axis.title = element_blank())
  
  ggsave(filename = paste0(COUNTRY,"-SOS-Seasonal.png"),
         plot = SOSmap,
         path = PlotDir,
         width= 180,
         height = 100,
         units = "mm",
         scale = 1,
         dpi = 600,
         type = "cairo",
         bg="white")
  
  ggsave(filename = paste0(COUNTRY,"-LGP-Seasonal.png"),
         plot = LGPmap,
         path = PlotDir,
         width= 90,
         height = 90,
         units = "mm",
         scale = 1,
         dpi = 600,
         type = "cairo",
         bg="white")
  
}



# Base plotting functions

add.map<-function(){plot(Map.Subset,add=T)}


Min<-min(Ystack[],na.rm=T)
Max<-max(Ystack[],na.rm=T)

plot(Ystack,
     type="classes",
     #col=rev( rainbow( 99, start=0,end=1 ) ), 
     breaks=seq(Min,Max,length.out=length(Min:Max)) ,
     col=pals::kovesi.cyclic_mygbm_30_95_c78_s25(n=36),
     fun=add.map)

plot(Ystack,
     type="continuous",
     col=pals::kovesi.cyclic_mygbm_30_95_c78_s25(n=36),
     fun=add.map)
