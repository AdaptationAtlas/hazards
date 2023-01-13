require(ggplot2)
require(data.table)
require(miceadds)
require(sp)
require(pbapply)
require(raster)
require(rworldmap)
require(rworldxtra)

# Where AI_Seasonal = T, I we end up with multiple seasons within a wet period for the LT_avg data. After the sequence merge logic is applied we should retain:
# 1) only the first sequence?
# 2) only the longest sequence?
# 3) all sequences, recoded to the same number (i.e. the season) with a field(s) that indicates the number and/or length of arid gaps.


# Import functions ####
source("R/chirps_sos/sos_functions.R")

# Set version ####
version<-1

DataDir<-"C:/Datasets"

# Set directories ####
CHIRPSraw_Dir<-paste0(DataDir,"/chirps_af/raw")
CHIRPS_Dekad_Dir<-paste0(DataDir,"/chirps_af/intermediate/array_countries_dekad")
SOS_Dir<-paste0(DataDir,"/atlas_SOS/intermediate/v",version)
list.dirs(SOS_Dir,recursive=F)
SOS_Dir<-list.dirs(SOS_Dir,recursive=F)[6]

# Load map of Africa ####
CHIRPSrast<-terra::rast(list.files(CHIRPSraw_Dir,full.names = T)[1])

AfricaMap<-rworldmap::getMap(resolution = "high")
AfricaMap<-AfricaMap[AfricaMap$REGION=="Africa"&!is.na(AfricaMap$REGION),]
AfricaMap<-sf::st_as_sf(AfricaMap)
AfricaMap<-terra::vect(AfricaMap)
AfricaMap<-terra::project(AfricaMap,CHIRPSrast)

Countries<-gsub("SOS_","",gsub(".RData","",list.files(SOS_Dir)))

# Wrangle SOS chunks ####
SOSFiles<-data.table(File=list.files(SOS_Dir,".RData",full.names = T),
                     Country=gsub("[.]RData","",tstrsplit(list.files(SOS_Dir,".RData",full.names = F),"_",keep=2)[[1]]))

AnalysisDir<-paste0(SOS_Dir,"/explore_results")
if(!dir.exists(AnalysisDir)){
  dir.create(AnalysisDir)
}

File<-paste0(AnalysisDir,"/SOSlist.RData")

# Join LTAvg data across countries
if(!file.exists(File)){
  SOSlist<-pblapply(1:nrow(SOSFiles),FUN=function(i){
    load(file=SOSFiles[i,File])
    list(LTAvg2=SOS_Data$LTAvg_SOS2,
         #Seasonal2=SOS_Data$Seasonal_SOS2,
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

NoSeasons<-41 #Total number of seasons in the LTavg dataset

MinLength<-1 # Minimum length of season allowed
Min.Rain<-0 # Minimum seasonal rainfall allowed
SeasonGapT<-3 # if the gap between seasons is <=SeasonGapT then it is considered as a split season (i.e. bimodal)

PlotDir1<-paste0(AnalysisDir,"/Plots/ML",MinLength,"-MR",Min.Rain,"-SG",SeasonGapT)

# The code will modify EOS, SOS and LGP according to SeasonGapT and attempt to tidy seasonal organization between S1 and S2

# LT avg ####
PlotDir<-paste0(PlotDir1,"/LT")
if(!dir.exists(PlotDir)){
  dir.create(PlotDir,recursive=T)
}


for(COUNTRY in Countries){
  print(COUNTRY)
  
  Map.Subset<-AfricaMap[AfricaMap$ADMIN==COUNTRY,]
  

   SOSData1<-data.table::copy(SOSlist[[COUNTRY]]$LTAVg_Summary)
   
   # NOTE - here seasons with insufficient length or rainfall are removed. This could be done after merging of bimodal seasons?
   SOSData1<-SOSData1[LGP>=MinLength & Tot.Rain>Min.Rain]
   
   if(grepl("AIsT",SOS_Dir)){
     
     select_seasons<-function(Seq,LGP,SOS){
       if(length(Seq)>2){
         ord<-order(LGP,decreasing = T)
         srt<-sort(LGP,decreasing = T)
         
         if(srt[2]==srt[3]){
           ord2<-ord[ord %in% c(2,3)]
           SOS<-SOS[ord2]
           SOSmax<-which(SOS==min(SOS))[1]
           ord<-c(1,ord2[SOSmax])
         }else{
           ord<-ord[ord %in% c(1,2)]
         }
         
       }else{
         return(rep(T,length(Seq)))
       }
     }
     
   }



   suppressWarnings(SOSData1[,Seq:=stringi::stri_replace_all_regex(Seq,
                                                                  pattern=unique(Seq[!is.na(Seq)]),
                                                                   replacement=1:length(unique(Seq[!is.na(Seq)]))),by=Index])
  
  
  # Load country raster
  CHIRPS_Ref<-terra::rast(paste0(CHIRPS_Dekad_Dir,"/",COUNTRY,".tif"))[[1]]
  CHIRPS_Ref<-data.table(terra::as.data.frame(CHIRPS_Ref,xy=T,cells=T))[,1:3]
  setnames(CHIRPS_Ref,"cell","Index")
  
  SOSData1<-merge(CHIRPS_Ref,SOSData1,by="Index",all.x=F)
  
  SOSData1<-dcast(SOSData1,Index+x+y+Seasons~Seq,value.var = c("SOS","EOS","LGP","Tot.Rain"))
  if(!"SOS_2" %in% colnames(SOSData1)){
    SOSData1[,c("SOS_2","EOS_2","LGP_2","Tot.Rain_2"):=NA]
  }
  
  # Calculate distance between seasons
  SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist12:=round(CicularDist(EOS_1,SOS_2,interval=36),0),by=Index]
  SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist21:=round(CicularDist(EOS_2,SOS_1,interval=36),0),by=Index]
  
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

# LT avg - Plot all countries ####
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

# LT avg - Plot all countries - SOS Maps ####
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

# LT avg - Plot all countries - LGP ####

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

# LT avg - Plot all countries - Seasonality Map ####

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

# Seasons only ####
Statistic<-"mean" # options are mean or mode

for(Statistic in c("mean","mode")){
  PlotDir<-paste0(PlotDir1,"/Seasonal-",Statistic)
  
  if(!dir.exists(PlotDir)){
    dir.create(PlotDir,recursive=T)
  }
  
  
  for(COUNTRY in Countries){
    print(paste0(COUNTRY,"-",Statistic))
    
    Map.Subset<-AfricaMap[AfricaMap$ADMIN==COUNTRY,]
    
      SOSData1<-data.table::copy(SOSlist[[COUNTRY]]$LTAvg2)
      setnames(SOSData1,c("Dekad.Season","Tot.Rain.mean"),c("Seq","Tot.Rain"))
      SOSData1[,Seasons:=length(unique(Seq)),by=Index]
      
      if(Statistic=="mean"){
        SOSData1[,SOS:=round(SOS.mean,0)][SOS==0,SOS:=36]
        SOSData1[,SOSdev:=SOS.mean.dev]
        SOSData1[,EOS:=round(EOS.mean,0)][EOS==0,EOS:=36]
        SOSData1[,EOSdev:=EOS.mean.dev]
        SOSData1[,LGP:=round(LGP.mean,0)]
      }
      
      if(Statistic=="mode"){
        SOSData1[,SOS:=round(SOS.mode,0)][SOS==0,SOS:=36]
        SOSData1[,SOSdev:=SOS.mode.dev]
        SOSData1[,EOS:=round(EOS.mode,0)][EOS==0,EOS:=36]
        SOSData1[,EOSdev:=EOS.mode.dev]
        SOSData1[,LGP:=round(LGP.mode,0)]
      }
  
    
    # NOTE - here seasons with insufficient length or rainfall are removed. This could be done after merging of bimodal seasons?
    SOSData1<-SOSData1[LGP>=MinLength & Tot.Rain>Min.Rain]
    
    # Load country raster
    CHIRPS_Ref<-terra::rast(paste0(CHIRPS_Dekad_Dir,"/",COUNTRY,".tif"))[[1]]
    CHIRPS_Ref<-data.table(terra::as.data.frame(CHIRPS_Ref,xy=T,cells=T))[,1:3]
    setnames(CHIRPS_Ref,"cell","Index")
    
    SOSData1<-merge(CHIRPS_Ref,SOSData1,by="Index",all.x=F)
    
    SOSData1<-dcast(SOSData1,Index+x+y+Seasons~Seq,value.var = c("Total.Seasons","SOS","SOSdev","EOS","EOSdev","LGP","LGP.sd","Tot.Rain","Tot.Rain.sd",
                                                               "Tot.ETo.mean","Tot.ETo.sd","Balance.mean","Balance.sd"))
    
    if(!"SOS_2" %in% colnames(SOSData1)){
      SOSData1[,c("Total.Seasons_2","SOS_2","SOSdev_2","EOS_2","EOSdev_2","LGP_2","LGP.sd_2","Tot.Rain_2","Tot.Rain.sd_2",
                  "Tot.ETo.mean_2","Tot.ETo.sd_2","Balance.mean_2","Balance.sd_2"):=NA]
    }
    
    # Calculate distance between seasons
    SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist12:=round(CicularDist(EOS_1,SOS_2,interval=36),0),by=Index]
    SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist21:=round(CicularDist(EOS_2,SOS_1,interval=36),0),by=Index]
    
    # Determine the order of seasons by the length of the distance between them, this is used to determine how bimodal
    # seasons are merged
    SOSData1[,MinDist:=paste(which(c(Dist12,Dist21)<=SeasonGapT),collapse="-"),by=Index][MinDist=="",MinDist:=NA]
    SOSData1[,MinDistN:=sum(c(Dist12,Dist21)<=SeasonGapT),by=Index]
    
    SOSData<-data.table::copy(SOSData1)
    
    # New columns S1, S2, etc. to record modified values
    SOSData[,S1:=SOS_1
            ][,S2:=SOS_2
              ][,SOSdev1:=SOSdev_1 
              ][,SOSdev2:=SOSdev_2 
              ][,E1:=EOS_1
              ][,EOSdev1:=EOSdev_1 
              ][,EOSdev2:=EOSdev_2 
                ][,E2:=EOS_2
                  ][,LGP1:=LGP_1
                    ][,LGP2:=LGP_2
                      ][,Tot.Rain1:=Tot.Rain_1
                        ][,Tot.Rain2:=Tot.Rain_2
                          ][,Tot.ETo1:=Tot.ETo.mean_1
                            ][,Tot.ETo2:=Tot.ETo.mean_2
                              ][,Balance1:=Balance.mean_1
                                ][,Balance2:=Balance.mean_2
                                  ][,Total.Seasons1:=Total.Seasons_1
                                    ][,Total.Seasons2:=Total.Seasons_2]
    
    # If both MinDist values are less than SeasonGapT then set end of season to dekad before SOS, unless S1==1 then set to 36
    SOSData[MinDistN==2,E1:=if(S1==1){36}else{S1-1},by=Index]  # Recalculate number of seasons using reclassified season cols
  
    # If there is one gap < MinDist and gap is between 1 & 2 then set end of season to E2
    SOSData[MinDistN==1 & MinDist==1,E1:=E2]
    SOSData[MinDistN==1 & MinDist==1,EOSdev1:=EOSdev2]
  
    # If there is one gap < MinDist and gap is between 2 & 1 then set start of season S1
    SOSData[MinDistN==1 & MinDist==2,S1:=S2]
    SOSData[MinDistN==1 & MinDist==2,SOSdev1:=SOSdev_2]
    
    # If there is one gap < MinDist and gap is between 2 & 1 then set start of season S1
    SOSData[MinDistN==1 & MinDist==2,Total.Seasons1:=max(Total.Seasons1,Total.Seasons2)]
  
    # For season we are combining sum LGP and then remove data relating to second season.
    SOSData[MinDistN %in% c(1,2),LGP1:=LGP1+LGP2
            ][MinDistN %in% c(1,2),Tot.Rain1:=Tot.Rain1+Tot.Rain2
              ][MinDistN %in% c(1,2),Tot.ETo1:=Tot.ETo1+Tot.ETo2
                ][MinDistN %in% c(1,2),Balance1:=Balance1+Balance2
                  ][MinDistN %in% c(1,2),c("LGP2","E2","S2","Tot.Rain2","Tot.ETo2","Balance2","EOSdev2","SOSdev2","Total.Seasons2"):=NA]
  
    SOSData[,Seasons:=sum(c(!is.na(S1),!is.na(S2)),na.rm=T),by=Index]
    
    # Create lists that correspond to dekads in each season
    if(SOSData[,!any(Seasons==2)]){ # If there is only 1 season and it is season 1
      SEASONS<-list(S1=OrderDekadSeq(unique(round(SOSData$S1))),S2=NULL)
      
      SOSData[,S1:=round(S1,0)][,S2:=NA][S1==0,S1:=36]
      
      SOSData[,Swettest:=S1][!is.na(S1),Seasons:=1]
    }else{
      if(SOSData[,!any(Seasons==1)]){ # If there is only 1 season and it is season 2
        SEASONS<-list(S1=OrderDekadSeq(unique(round(SOSData$S2))),S2=NULL)
        
        SOSData[,S1:=round(S1,0)][,S2:=NA][S1==0,S1:=36]
        
        SOSData[,Swettest:=S1][!is.na(S1),Seasons:=1]
      }else{ # Multiple seasons
        # Use multimode function to determine seasons from the distribution of SOS values (from all cells and both seasons combined)
        # This approach doesn't work well for some countries for which their seasons are specified manually
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
        
        # Reorganize seasons - create columns to state which season a cell x season belongs to
        SOSData[,X1:=as.numeric(NA)
        ][,X2:=as.numeric(NA)
        ][S1 %in% SEASONS$S1,X1:=S1
        ][S1 %in% SEASONS$S2,X2:=S1
        ][S2 %in% SEASONS$S1,X1:=S2
        ][S2 %in% SEASONS$S2,X2:=S2
        ][S1 %in% SEASONS$S1 & S2 %in% SEASONS$S1,X1:=S1 # Where both S1 and S2 are in the same season then assign one to S1 and the other to S2
        ][S1 %in% SEASONS$S1 & S2 %in% SEASONS$S1,X2:=S2 
        ][S1 %in% SEASONS$S2 & S2 %in% SEASONS$S2,X1:=S1
        ][S1 %in% SEASONS$S2 & S2 %in% SEASONS$S2,X2:=S2]
  
        # Assign LGP to reorganised seasons
        SOSData[X1==S1,Y1:=LGP1
                ][X1==S2,Y1:=LGP2
                  ][X2==S1,Y2:=LGP1
                    ][X2==S2,Y2:=LGP2
                      ][,LGP1:=Y1
                        ][,LGP2:=Y2
                          ][,Y1:=NULL
                            ][,Y2:=NULL]
        
        # Assign SOSdev to reorganised seasons
        SOSData[X1==S1,Y1:=SOSdev1
        ][X1==S2,Y1:=SOSdev2
        ][X2==S1,Y2:=SOSdev1
        ][X2==S2,Y2:=SOSdev2
        ][,SOSdev1:=Y1
        ][,SOSdev2:=Y2
        ][,Y1:=NULL
        ][,Y2:=NULL]
        
        # Assign EOSdev to reorganised seasons
        SOSData[X1==S1,Y1:=EOSdev1
        ][X1==S2,Y1:=EOSdev2
        ][X2==S1,Y2:=EOSdev1
        ][X2==S2,Y2:=EOSdev2
        ][,EOSdev1:=Y1
        ][,EOSdev2:=Y2
        ][,Y1:=NULL
        ][,Y2:=NULL]
        
        # Assign Tot.Rain to reorganised seasons
        SOSData[X1==S1,Y1:=Tot.Rain1
        ][X1==S2,Y1:=Tot.Rain2
        ][X2==S1,Y2:=Tot.Rain1
        ][X2==S2,Y2:=Tot.Rain2
        ][,Tot.Rain1:=Y1
        ][,Tot.Rain2:=Y2
        ][,Y1:=NULL
        ][,Y2:=NULL]
        
        # Assign Tot.ETo to reorganised seasons
        SOSData[X1==S1,Y1:=Tot.ETo1
        ][X1==S2,Y1:=Tot.ETo2
        ][X2==S1,Y2:=Tot.ETo1
        ][X2==S2,Y2:=Tot.ETo2
        ][,Tot.ETo1:=Y1
        ][,Tot.ETo2:=Y2
        ][,Y1:=NULL
        ][,Y2:=NULL]
        
        # Assign Balance to reorganised seasons
        SOSData[X1==S1,Y1:=Balance1
        ][X1==S2,Y1:=Balance2
        ][X2==S1,Y2:=Balance1
        ][X2==S2,Y2:=Balance2
        ][,Balance1:=Y1
        ][,Balance2:=Y2
        ][,Y1:=NULL
        ][,Y2:=NULL]
        
        # Assign Total.Seasons to reorganised seasons
        SOSData[X1==S1,Y1:=Total.Seasons1
        ][X1==S2,Y1:=Total.Seasons2
        ][X2==S1,Y2:=Total.Seasons1
        ][X2==S2,Y2:=Total.Seasons2
        ][,Total.Seasons1:=Y1
        ][,Total.Seasons2:=Y2
        ][,Y1:=NULL
        ][,Y2:=NULL]
        
        # Assign EOS to reorganised seasons
        SOSData[X1==S1,Y1:=E1
        ][X1==S2,Y1:=E2
        ][X2==S1,Y2:=E1
        ][X2==S2,Y2:=E2
        ][,E1:=Y1
        ][,E2:=Y2
        ][,Y1:=NULL
        ][,Y2:=NULL]
        
        SOSData[,S1:=X1][,S2:=X2][,X1:=NULL][,X2:=NULL]
  
        # Assign wettest season
        SOSData[,Swettest:=S1
        ][!is.na(S1) & Tot.Rain1>Tot.Rain2,Swettest:=S2
        ][is.na(Swettest),Swettest:=S2]
        
      }}
    
    SOSData[,Seasons:=sum(c(!is.na(S1),!is.na(S2))),by=Index]
    
    # Manually tidy up countries which are largely unimodal but have a few pixels that split the seasons.
    # Note: Check that correct LGP field is being used here
    if(COUNTRY %in% c("Nigeria","Burundi","Madagascar","Mozambique")){
      SOSData[is.na(S1) & !is.na(S2),E1:=E2
      ][is.na(S1) & !is.na(S2),E2:=NA
      ][is.na(S1) & !is.na(S2),LGP1 :=LGP2
      ][is.na(S1) & !is.na(S2),LGP2:=NA
      ][is.na(S1) & !is.na(S2),SOSdev1:=SOSdev2
      ][is.na(S1) & !is.na(S2),SOSdev2:=NA
      ][is.na(S1) & !is.na(S2),EOSdev1:=EOSdev2
      ][is.na(S1) & !is.na(S2),EOSdev2:=NA
      ][is.na(S1) & !is.na(S2),Tot.Rain1:=Tot.Rain2
      ][is.na(S1) & !is.na(S2),Tot.Rain2:=NA
      ][is.na(S1) & !is.na(S2),Tot.ETo1:=Tot.ETo2
      ][is.na(S1) & !is.na(S2),Tot.ETo2:=NA
      ][is.na(S1) & !is.na(S2),Balance1:=Balance2
      ][is.na(S1) & !is.na(S2),Balance2:=NA
      ][is.na(S1) & !is.na(S2),Total.Seasons1:=Total.Seasons2
      ][is.na(S1) & !is.na(S2),Total.Seasons2:=NA
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
    
    # TO DO: 2-factor pallete for deviance x magnitude - could be adjusting alpha or a biscale plot ####
    
    # SOS Map Data ####
    S1<-round(terra::rasterize(X,BaseRaster,field="S1"),0)
    if(SOSData[!is.na(S2),.N==0]){
      S2<-S1
      S2[]<-NA
    }else{
      S2<-round(terra::rasterize(X,BaseRaster,field="S2"),0)
    }
    
    Swettest<-round(terra::rasterize(X,BaseRaster,field="Swettest"),0)
    
    SOS<-c(S1,S2,Swettest)
    
    names(SOS)<-c("SOS_1",
                  "SOS_2",
                  "SOS_Wettest")
    
    SOS2<-raster::stack(SOS)
    names(SOS2) <- names(SOS)
    
    coords <- xyFromCell(SOS2, seq_len(ncell(SOS2)))
    
    # SOS Deviance Map Data ####
    SOSdev1<-terra::rasterize(X,BaseRaster,field="SOSdev1")
    if(SOSData[!is.na(SOSdev2),.N==0]){
      SOSdev2<-SOSdev1
      SOSdev2[]<-NA
    }else{
      SOSdev2<-terra::rasterize(X,BaseRaster,field="SOSdev2")
    }
    
    SOSdev<-c(SOSdev1,SOSdev2)
    
    names(SOSdev)<-c("SOSdev1","SOSdev2")
    
    SOSdev2<-raster::stack(SOSdev)
    names(SOSdev2) <- names(SOSdev)
    
    # LGP Map Data ####
    LGP_1<-round(terra::rasterize(X,BaseRaster,field="LGP1"),0)
    
    if(!2 %in% SOSData$Seasons){
      LGP_2<-LGP_1
      LGP_2[]<-NA
    }else{
      LGP_2<-round(terra::rasterize(X,BaseRaster,field="LGP2"),0)
    }
    
    LGP<-c(LGP_1,LGP_2,sum(c(LGP_1,LGP_2),na.rm = T))
    
    names(LGP)<-c("LGP_1",
                  "LGP_2",
                  "LGP_sum")
    
    LGP2<-raster::stack(LGP)
    names(LGP2) <- names(LGP)
    
    
    # Total Rain Map Data ####
    TotRain1<-round(terra::rasterize(X,BaseRaster,field="Tot.Rain1"),0)
    
    if(!2 %in% SOSData$Seasons){
      TotRain2<-TotRain1
      TotRain2[]<-NA
    }else{
      TotRain2<-round(terra::rasterize(X,BaseRaster,field="Tot.Rain2"),0)
    }
    
    TotRain<-c(TotRain1,TotRain2,sum(c(TotRain1,TotRain2),na.rm = T))
    
    names(TotRain)<-c("TotRain1",
                  "TotRain2",
                  "TotRain_sum")
    
    TotRain2<-raster::stack(TotRain)
    names(TotRain2) <- names(TotRain)
    
    # Balance Map Data ####
    Balance1 <-round(terra::rasterize(X,BaseRaster,field="Balance1"),0)
    
    if(!2 %in% SOSData$Seasons){
      Balance2<-Balance1 
      Balance2[]<-NA
    }else{
      Balance2<-round(terra::rasterize(X,BaseRaster,field="Balance2"),0)
    }
    
    Balance<-c(Balance1 ,Balance2,sum(c(Balance1 ,Balance2),na.rm = T))
    
    names(Balance)<-c("Balance1 ",
                      "Balance2",
                      "Balance2_sum")
    
    Balance2<-raster::stack(Balance)
    names(Balance2) <- names(Balance)
    
    # Total Seasons Map Data ####
    Total.Seasons1 <-terra::rasterize(X,BaseRaster,field="Total.Seasons1")
    
    if(!2 %in% SOSData$Seasons){
      Total.Seasons2<-Total.Seasons1 
      Total.Seasons2[]<-NA
    }else{
      Total.Seasons2<-terra::rasterize(X,BaseRaster,field="Total.Seasons2")
    }
    
    Total.Seasons<-c(Total.Seasons1 ,Total.Seasons2)
    Total.Seasons<-Total.Seasons/NoSeasons
    
    names(Total.Seasons)<-c("Prop.Seasons1","Prop.Seasons2")                  
    
    Total.Seasons2<-raster::stack(Total.Seasons)
    names(Total.Seasons2) <- names(Total.Seasons)
    
    # Plot SOS ####
    Data <- raster::stack(as.data.frame(raster::getValues(SOS2)))
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
    
    # Plot SOS deviance ####
    Data <- raster::stack(as.data.frame(raster::getValues(SOSdev2)))
    
    Data <- cbind(coords, Data)
  
    SOSdevmap<-ggplot(Data) + 
      geom_tile(aes(x, y, fill = values)) +
      facet_wrap(~ ind,ncol=2) +
      viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
      theme_bw()+
      labs(fill="SOS deviance (dekads)")+
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            legend.key.width=unit(0.1,"npc"))+
      coord_equal()
    
    # Plot LGP ####
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
            axis.title = element_blank(),
            legend.key.width=unit(0.1,"npc"))+
      coord_equal()
    
    # Plot Rain ####
    Data <- raster::stack(as.data.frame(raster::getValues(TotRain2)))
    
    Data <- cbind(coords, Data)
    
    Rainmap<-ggplot(Data) + 
      geom_tile(aes(x, y, fill = values)) +
      facet_wrap(~ ind,ncol=3) +
      viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
      theme_bw()+
      labs(fill="Total Rainfall (mm)")+
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            legend.key.width=unit(0.1,"npc"))+
      coord_equal()
    
    # Plot Balance ####
    Data <- raster::stack(as.data.frame(raster::getValues(Balance2)))
    
    Data <- cbind(coords, Data)
    
    Balancemap<-ggplot(Data) + 
      geom_tile(aes(x, y, fill = values)) +
      facet_wrap(~ ind,ncol=3) +
      viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
      theme_bw()+
      labs(fill="Total Rainfall - Total ETo (mm)")+
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            legend.key.width=unit(0.1,"npc"))+
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
    
    Seasonalitymap<-ggplot(Data) + 
      geom_tile(aes(x, y, fill = System)) +
      viridis::scale_fill_viridis(option="turbo",discrete = T,na.value="transparent",na.translate = F)+
      theme_bw()+
      labs(fill="System")+
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            legend.direction = "horizontal")+
      guides(fill=guide_legend(ncol=1))+
      coord_equal()
    
    # Number Seasons Map ####
    Yseasons<-terra::rasterize(X,BaseRaster,field="Seasons")
    names(Yseasons)<-"Seasons"
    coords <- xyFromCell(Yseasons, seq_len(ncell(Yseasons)))
    Data <-values(Yseasons)
    names(Yseasons) <- names(Yseasons)
    
    Data <- data.frame(cbind(coords, Seasons=Data))
    
    Noseasons<-ggplot(Data) + 
      geom_tile(aes(x, y, fill = as.factor(Seasons))) +
      viridis::scale_fill_viridis(option="viridis",discrete = T,na.value="transparent",na.translate = F,direction=-1)+
      theme_bw()+
      labs(fill="Seasons")+
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            legend.direction = "horizontal")+
      guides(fill=guide_legend(ncol=2))+
      coord_equal()
    # Failed Seasons Map ####
    Data <- raster::stack(as.data.frame(100*(1-raster::getValues(Total.Seasons2))))
    
    Data <- cbind(coords, Data)
    
    Failedmap<-ggplot(Data) + 
      geom_tile(aes(x, y, fill = values)) +
      facet_wrap(~ ind,ncol=3) +
      viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
      theme_bw()+
      labs(fill="% Failed Seasons")+
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            legend.key.width=unit(0.1,"npc"))+
      coord_equal()
    
    # Combine plots and save ####
    
  
    g1<-gridExtra::arrangeGrob(grobs=list(SOSmap,LGPmap,Rainmap),layout_matrix=matrix(rep(1:3,3),nrow = 3))
    g2<-gridExtra::arrangeGrob(grobs=list(SOSdevmap,Failedmap,Balancemap),layout_matrix=matrix(rep(1:3,3),nrow = 3))
    g3<-gridExtra::arrangeGrob(grobs=list(Noseasons,Seasonalitymap),layout_matrix=matrix(c(1,1,2,2,2),ncol = 1))
      
    
    gcomb<-gridExtra::arrangeGrob(grobs=list(g3,g1,g2),layout_matrix=matrix(c(1,2,2,2,3,3,3),nrow = 1))
    
    ggsave(filename = paste0(COUNTRY,".png"),
           plot = gcomb,
           path = PlotDir,
           width= 200,
           height = 200*0.6,
           units = "mm",
           scale = 2.5,
           device=png,
           dpi = 600,
           bg="white")
    
  }
  
  # Seasonal - Plot all countries ####
  Files<-list.files(paste0(PlotDir,"/Data"),".RData")
  
  Data<-miceadds::load.Rdata2(filename=Files[1],path=paste0(PlotDir,"/Data"))
  
  SOSData<-data.table::rbindlist(lapply(Files,FUN=function(X){miceadds::load.Rdata2(filename=X,path=paste0(PlotDir,"/Data"))}),use.names = T)
  SOSData[,Seasons:=sum(!is.na(SOS_1),!is.na(SOS_2)),by=list(x,y,Index)]
  
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
  
  
  # Seasonal - SOS Map Data ####
  S1<-round(terra::rasterize(X,BaseRaster,field="SOS_1"),0)
  if(SOSData[!is.na(S2),.N==0]){
    S2<-S1
    S2[]<-NA
  }else{
    S2<-round(terra::rasterize(X,BaseRaster,field="SOS_2"),0)
  }
  
  Swettest<-round(terra::rasterize(X,BaseRaster,field="Swettest"),0)
  
  SOS<-c(S1,S2,Swettest)
  
  names(SOS)<-c("SOS_1",
                "SOS_2",
                "SOS_Wettest")
  
  SOS2<-raster::stack(SOS)
  names(SOS2) <- names(SOS)
  
  coords <- xyFromCell(SOS2, seq_len(ncell(SOS2)))
  
  # Seasonal - SOS Deviance Map Data ####
  SOSdev1<-terra::rasterize(X,BaseRaster,field="SOSdev_1")
  if(SOSData[!is.na(SOSdev2),.N==0]){
    SOSdev2<-SOSdev1
    SOSdev2[]<-NA
  }else{
    SOSdev2<-terra::rasterize(X,BaseRaster,field="SOSdev_2")
  }
  
  SOSdev<-c(SOSdev1,SOSdev2)
  
  names(SOSdev)<-c("SOSdev1","SOSdev2")
  
  SOSdev2<-raster::stack(SOSdev)
  names(SOSdev2) <- names(SOSdev)
  
  # Seasonal - LGP Map Data ####
  LGP_1<-round(terra::rasterize(X,BaseRaster,field="LGP_1"),0)
  LGP_2<-round(terra::rasterize(X,BaseRaster,field="LGP_2"),0)
  
  LGP<-c(LGP_1,LGP_2,sum(c(LGP_1,LGP_2),na.rm = T))
  
  names(LGP)<-c("LGP_1",
                "LGP_2",
                "LGP_sum")
  
  LGP2<-raster::stack(LGP)
  names(LGP2) <- names(LGP)
  
  # Seasonal - Total Rain Map Data ####
  TotRain1<-round(terra::rasterize(X,BaseRaster,field="Tot.Rain_1"),0)
  TotRain2<-round(terra::rasterize(X,BaseRaster,field="Tot.Rain_2"),0)
  TotRain<-c(TotRain1,TotRain2,sum(c(TotRain1,TotRain2),na.rm = T))
  
  names(TotRain)<-c("TotRain1",
                    "TotRain2",
                    "TotRain_sum")
  
  TotRain2<-raster::stack(TotRain)
  names(TotRain2) <- names(TotRain)
  
  # Seasonal - Balance Map Data ####
  Balance1 <-round(terra::rasterize(X,BaseRaster,field="Balance.mean_1"),0)
  Balance2<-round(terra::rasterize(X,BaseRaster,field="Balance.mean_2"),0)
  Balance<-c(Balance1 ,Balance2,sum(c(Balance1 ,Balance2),na.rm = T))
  
  names(Balance)<-c("Balance1 ",
                    "Balance2",
                    "Balance2_sum")
  
  Balance2<-raster::stack(Balance)
  names(Balance2) <- names(Balance)
  
  # Seasonal - Total Seasons Map Data ####
  Total.Seasons1 <-terra::rasterize(X,BaseRaster,field="Total.Seasons_1")
  Total.Seasons2<-terra::rasterize(X,BaseRaster,field="Total.Seasons_2")
  Total.Seasons<-c(Total.Seasons1 ,Total.Seasons2)
  Total.Seasons<-Total.Seasons/NoSeasons
  
  names(Total.Seasons)<-c("Prop.Seasons1","Prop.Seasons2")                  
  
  Total.Seasons2<-raster::stack(Total.Seasons)
  names(Total.Seasons2) <- names(Total.Seasons)
  
  # Seasonal - Plot SOS ####
  Data <- raster::stack(as.data.frame(raster::getValues(SOS2)))
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
         height = 200*0.6,
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Plot SOS deviance ####
  Data <- raster::stack(as.data.frame(raster::getValues(SOSdev2)))
  
  Data <- cbind(coords, Data)
  
  SOSdevmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=2) +
    viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
    theme_bw()+
    labs(fill="SOS deviance (dekads)")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.key.width=unit(0.1,"npc"))+
    coord_equal()
  
  ggsave(filename = "Africa - SOSdev.png",
         plot = SOSdevmap,
         path = PlotDir,
         width= 200,
         height = 200*0.6,      
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Plot LGP ####
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
          axis.title = element_blank(),
          legend.key.width=unit(0.1,"npc"))+
    coord_equal()
  
  ggsave(filename = "Africa - LGP.png",
         plot = LGPmap,
         path = PlotDir,
         width= 200,
         height = 200*0.6,      
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Plot Rain ####
  Data <- raster::stack(as.data.frame(raster::getValues(TotRain2)))
  
  Data <- cbind(coords, Data)
  
  Rainmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=3) +
    viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
    theme_bw()+
    labs(fill="Total Rainfall (mm)")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.key.width=unit(0.1,"npc"))+
    coord_equal()
  
  ggsave(filename = "Africa - Rain.png",
         plot = Rainmap,
         path = PlotDir,
         width= 200,
         height = 200*0.6,      
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Plot Balance ####
  Data <- raster::stack(as.data.frame(raster::getValues(Balance2)))
  
  Data <- cbind(coords, Data)
  
  Balancemap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=3) +
    viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
    theme_bw()+
    labs(fill="Total Rainfall - Total ETo (mm)")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.key.width=unit(0.1,"npc"))+
    coord_equal()
  
  ggsave(filename = "Africa - Balance.png",
         plot = Balancemap,
         path = PlotDir,
         width= 200,
         height = 200*0.6,      
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Plot Seasonality ####
  
  Ysystem<-terra::rasterize(X,BaseRaster,field="System")
  names(Ysystem)<-"System"
  coords <- xyFromCell(Ysystem, seq_len(ncell(Ysystem)))
  Data <-values(Ysystem)
  names(Ysystem) <- names(Ysystem)
  
  Levels<-data.table(value=sort(unique(Ysystem[!is.na(Ysystem)])),description=levels(Ysystem)[[1]])
  
  Data <- data.frame(cbind(coords, Data))
  Data$System<-Levels$description.category[match(Data$System,Levels$value)]
  
  Seasonalitymap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = System)) +
    viridis::scale_fill_viridis(option="turbo",discrete = T,na.value="transparent",na.translate = F)+
    theme_bw()+
    labs(fill="System")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.direction = "horizontal")+
    guides(fill=guide_legend(ncol=1))+
    coord_equal()
  
  ggsave(filename = "Africa - Seasonality.png",
         plot = Seasonalitymap,
         path = PlotDir,
         width= 200,
         height = 200,      
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Plot Number Seasons ####
  Yseasons<-terra::rasterize(X,BaseRaster,field="Seasons")
  names(Yseasons)<-"Seasons"
  coords <- xyFromCell(Yseasons, seq_len(ncell(Yseasons)))
  Data <-values(Yseasons)
  names(Yseasons) <- names(Yseasons)
  
  Data <- data.frame(cbind(coords, Seasons=Data))
  
  Noseasons<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = as.factor(Seasons))) +
    viridis::scale_fill_viridis(option="viridis",discrete = T,na.value="transparent",na.translate = F,direction=-1)+
    theme_bw()+
    labs(fill="Seasons")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.direction = "horizontal")+
    guides(fill=guide_legend(ncol=2))+
    coord_equal()
  
  ggsave(filename = "Africa - NoSeasons.png",
         plot = Noseasons,
         path = PlotDir,
         width= 200,
         height = 200,      
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Plot Failed Seasons ####
  Data <- raster::stack(as.data.frame(100*(1-raster::getValues(Total.Seasons2))))
  
  Data <- cbind(coords, Data)
  
  Failedmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=3) +
    viridis::scale_fill_viridis(option="viridis",discrete = F,direction = -1,na.value = "transparent") +
    theme_bw()+
    labs(fill="% Failed Seasons")+
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          legend.key.width=unit(0.1,"npc"))+
    coord_equal()
  
  ggsave(filename = "Africa - FailedSeasons.png",
         plot = Failedmap,
         path = PlotDir,
         width= 200,
         height = 200*0.6,      
         units = "mm",
         scale = 1.2,
         dpi = 600,
         device = png,
         bg="white")
  
  # Seasonal - Combine plots and save ####
  
  
  g1<-gridExtra::arrangeGrob(grobs=list(SOSmap,LGPmap,Rainmap),layout_matrix=matrix(rep(1:3,3),nrow = 3))
  g2<-gridExtra::arrangeGrob(grobs=list(SOSdevmap,Failedmap,Balancemap),layout_matrix=matrix(rep(1:3,3),nrow = 3))
  g3<-gridExtra::arrangeGrob(grobs=list(Noseasons,Seasonalitymap),layout_matrix=matrix(c(1,1,2,2,2),ncol = 1))
  
  
  gcomb<-gridExtra::arrangeGrob(grobs=list(g3,g1,g2),layout_matrix=matrix(c(1,2,2,2,3,3,3),nrow = 1))
  
  ggsave(filename = "Africa - Combined.png",
         plot = gcomb,
         path = PlotDir,
         width= 200,
         height = 200*0.6,
         units = "mm",
         scale = 2.5,
         device=png,
         dpi = 600,
         bg="white")
}