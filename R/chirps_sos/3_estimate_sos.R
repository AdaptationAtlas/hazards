require(lubridate)
require(data.table)
require(miceadds)
require(future)
require(terra)
require(circular)
require(sp)

options(scipen=999)

#Cores<-parallel::detectCores()
Cores<-10

# Data and save locations  ####
#DataDir<-"/home/jovyan/common_data"
DataDir<-"C:/Datasets"

# Import Functions ####
source("R/chirps_sos/sos_functions.R")

# Version
version<-1

# Set directories
CHIRPSraw_Dir<-paste0(DataDir,"/chirps_af/raw")
Hobbins_Dir<-paste0(DataDir,"/hobbins_ref_et/intermediate/array_countries")
CHIRPS_Dekad_Dir<-paste0(DataDir,"/chirps_af/intermediate/array_countries_dekad")
SOS_Dir<-paste0(DataDir,"/atlas_SOS/intermediate/v",version)

if(!dir.exists(SOS_Dir)){
  dir.create(SOS_Dir,recursive = T)
}

if(!dir.exists(CHIRPS_Dekad_Dir)){
  dir.create(CHIRPS_Dekad_Dir)
}

# Set Parameters ####
Debug<-F
rmNAs<-T

# Set Minimum rainfall required for second season
Min.Rain<-0

# Set separation distance in months of first and second season
# SeasonSeparation<-2 # This is not currently implemented?

# For a 150 day season planting starts 7.5 dekads before mid-point of three month wet period (e.g. dekad 4.5 - 7.5 = -3)
# We therefore pad the rainy month periods by 3 dekads in each direction
PadBack<-3
PadForward<-3

Do.SeqMerge.LT<-T # Merge sequences for climatological classification
  
D1.mm<-25 # Amount of rainfall needed in dekad n for SOS to be TRUE 
D2.mm<-20 # Amount of rainfall needed in dekad n+1 + dekad n+2 for SOS to be TRUE in dekad n 
D2.len<-2 # Number of dekads following dekad n for which rain is summed and assessed against the threshold presented in D2.mm 
AI.t<-0.5 # Aridity index threshold that dekadal rainfall must fall below to indicate EOS
Do.SeqMerge<-T # Merge sequences that are close together and remove false starts
MaxGap<-1 # An integer value describing the maximum gap (number of NA values) allowed between non-NA values before the sequence breaks.
MinStartLen<-1 # An integer value describing the minimum length of first sequence block, if the first sequence is too short and too separated from the next sequence it is removed. This is to remove false starts.
MaxStartSep<-1 # an integer value describing the maximum separation of the first sequence block, if the first sequence is too short and too separated from the next sequence it is removed. This is to remove false starts.
ClipAI<-F # logical `T/F`, if `T` then `AI` values corresponding to the last non-NA value in `Seq` are all set to `F` and the sequence is halted by the last `F` value of AI.
Season2.Prop<-0.25 # Proportions of seasons the second or third seasons need to be present so that they are included in outputs.
MinLength<-1 # Minimum length of second or third growing season in dekads
AI_Seasonal<-F # Calculate aridity index of each dekad on an annual basis (T) or using the long-term average across all years (F)
RollBack<-T

S1.AI<-T # Set AI to T when for the first value of sequences of RAIN == T


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

FolderName<-paste(c(
  "S2mm",Min.Rain,
  "-Pad",PadBack,"x",PadForward,
  "-D1mm",D1.mm,
  "-D2mm",D2.mm,
  "-D2l",D2.len,
  "-AIt",AI.t, 
  "-SqM",substr(Do.SeqMerge,1,1),
  "-SqMLT",substr(Do.SeqMerge.LT,1,1),
  "-MxGp",MaxGap,
  "-MiSL",MinStartLen,
  "-MaSL",MaxStartSep,
  "-ClAI",substr(ClipAI,1,1),
  "-S2Pr",Season2.Prop,
  "-S2l",MinLength,
  "-AIs",substr(AI_Seasonal,1,1),
  "-RB",substr(RollBack,1,1),
  "-S1AI",S1.AI
),collapse = "")

SaveDir<-paste0(SOS_Dir,"/",FolderName)

if(!dir.exists(SaveDir)){
  dir.create(SaveDir,recursive=T)
}

Overwrite<-T # Overwrite previous saved analysis?

# Loop SOS analysis over chunks ####
for(COUNTRY in Countries){
  FILE<-paste0(SaveDir,"/SOS_",COUNTRY,".RData")
  if(!file.exists(FILE) | Overwrite==T){
  print(paste("Loading Country",COUNTRY," - ", match(COUNTRY,Countries),"/",length(Countries)))
  
  # 1) Load & Prepare Data #####
  
  # Load CHIRPS
  CHIRPSrast<-terra::rast(paste0(CHIRPS_Dekad_Dir,"/",COUNTRY,".tif"))

  # Load ETo
  HOBBINSrast<-terra::rast(paste0(Hobbins_Dir,"/",COUNTRY,".tif"))
  # Reformat ETo names to match CHIRPS
  names(HOBBINSrast)<-gsub("pet_","",names(HOBBINSrast))
  names(HOBBINSrast)<-paste0(substr(names(HOBBINSrast),1,4),".",as.numeric(substr(names(HOBBINSrast),5,6)))
  
  # Make sure timeframe of ETo and CHIRPS are aligned
  HOBBINSrast<-HOBBINSrast[[names(HOBBINSrast) %in% names(CHIRPSrast)]]
  CHIRPSrast<-CHIRPSrast[[names(CHIRPSrast) %in% names(HOBBINSrast)]]
  
  # 1.1) Convert raster stacks to tables ######
  
  # Rows are cells and columns layers (i.e. dekads) 
  CHIRPS<-data.table(terra::as.data.frame(CHIRPSrast,xy=T,cells=T))
  HOBBINS<-data.table(terra::as.data.frame(HOBBINSrast,xy=T,cells=T))
  
  # reformat to make dataset long each row is a cell x dekad (lyr) value
  CLIM<-melt(CHIRPS,id.vars=c("cell","x","y"),value.name = "Rain")
  HOBBINS<-melt(HOBBINS,id.vars=c("cell","x","y"),value.name = "ETo")
  
  # Check cells indices match
  print(paste0("All cells indices between CHIRPS and HOBBINS match? = ",all(CLIM$cell == HOBBINS$cell)))

  # Add ETo to CHIRPS table
  CLIM[,ETo:=HOBBINS[,ETo]
       ][,Dekad:=as.integer(unlist(data.table::tstrsplit(variable,"[.]",keep=2))),by=variable
         ][,Year:=as.integer(unlist(data.table::tstrsplit(variable,"[.]",keep=1))),by=variable
           ][,Month:=as.integer(dkd2mth(Dekad)),by=Dekad
             ][,variable:=NULL
               ][,AI:=round(Rain/ETo,2)   # Calculate aridity index (AI)
                 ][,Rain:=round(Rain,2)
                   ][,ETo:=round(ETo,2)]
  
  setnames(CLIM,"cell","Index")
  
  # Remove NA values & order
  CLIM<-CLIM[!(ETo<0 | Rain<0)][order(Index,Year,Dekad)]
  
  # Check for incomplete data
  Incomplete<-CLIM[,.N,by=Index][N<max(N)]
  if(nrow(Incomplete)>0){
    print("These cells have incomplete data:")
    print(Incomplete)
  }
  
  CLIM<-CLIM[!Index %in% Incomplete$Index]
  
  # Run below to check for any missing values
  if(F){
  xy<-unique(CLIM[,list(x,y)])
  terra::plot(terra::vect(xy,geom=c("x","y")))
  }
  
  # Tidy up
  rm(HOBBINS,CHIRPSrast,HOBBINSrast,Incomplete)
  gc()
  

  format(object.size(CLIM),units="Gb")
  
  # 2) Calculate seasonality using long-term average data only #####
  print(paste("Calculating LT Avg:",COUNTRY,"-", match(COUNTRY,Countries),"/",length(Countries)))
  
  # The below can made parallel
  CLIM.LT<-data.table::copy(CLIM)[,Rain.sum2:=slide_apply(data=Rain,window=D2.len+1,step=1,fun=sum),by=list(Index)
                                  ][,Rain.sum9:=slide_apply2(Rain,window=9,step=1,fun=sum),by=list(Index) # Rainfall over 9 dekads (3 months)
                                      ][,list(ETo=mean(ETo,na.rm=T),
                                              Rain=mean(Rain,na.rm=T),
                                              Rain.sum2=mean(Rain.sum2,na.rm=T),
                                              Rain.sum9=mean(Rain.sum9,na.rm=T)),by=list(Index,Dekad)
                                        ][,AI:=round(Rain/ETo,2)
                                         ][,AImet:=AI>=AI.t
                                           ][,SOSmet:=Rain.sum2>=D2.mm & Rain>=D1.mm
                                            ][,Seq2:=SOS_RSeasonLT(RAIN=SOSmet,AI=AImet,S1.AI=S1.AI),by=list(Index)]
  
  
  # Merge sequences over year end boundary? (should always TRUE, I think option was to quality control LTSeqMerge function)
  if(Do.SeqMerge.LT){
    CLIM.LT[,Seq:=LTSeqMerge(Seq2),by=Index]
  }else{
    CLIM.LT[,Seq:=Seq2]
  }
  
  if(CLIM.LT[,!all(is.na(Seq))]){
  CLIM.LT[!is.na(Seq),SOS:=OrderDekadSeq(Dekad)[1],by=list(Index,Seq)
          ][!is.na(Seq),EOS:= tail(OrderDekadSeq(Dekad), n=1),by=list(Index,Seq)
            ][!is.na(Seq),LGP:=.N,by=list(Index,Seq)
              ][!is.na(Seq),Tot.Rain:=sum(Rain),by=list(Index,Seq)
                ][!is.na(Seq),Tot.ETo:=sum(ETo),by=list(Index,Seq)
                  ][,AImet_sum:=sum(AImet),by=Index]
  
  # Where AI is met in all 36 dekads, choose the first instance where SOSmet is T in the 3 wettest months
  CLIM.LT[AImet_sum==36,SOS:=round(median(Dekad[Rain.sum9==max(Rain.sum9)])),by=Index
          ][AImet_sum==36,EOS:=if(SOS[1]==1){36}else{SOS[1]-1},by=Index]
  }else{
    CLIM.LT[,SOS:=as.numeric(NA)
            ][,EOS:=as.numeric(NA) 
              ][,LGP:=as.numeric(NA)
                ][,Tot.Rain:=as.numeric(NA)
                  ][,Tot.ETo:=as.numeric(NA)
                    ][,AImet_sum:=as.numeric(NA)]
  }
  
  CLIM.LT.Summary<-unique(CLIM.LT[!is.na(Seq),list(Index,Seq,SOS,EOS,LGP,Tot.Rain,Tot.ETo,AImet_sum)])
  CLIM.LT.Summary[,Balance:=Tot.Rain/Tot.ETo][,Seasons:=.N,by=Index]
  
  # 3) Find wet months #####
  # This section could perhaps be moved to the SOS_Fun function?
  print(paste("Calculating Wet Months:",COUNTRY,"-", match(COUNTRY,Countries),"/",length(Countries)))
  
  CLIM.MONTH<-CLIM[,list(Rain.M.Sum=sum(Rain)),by=list(Index,Year,Month) # Sum rainfall by month within year and site
                   ][,Rain.M.Sum3:=slide_apply2(Rain.M.Sum,window=3,step=1,fun=sum),by=list(Index) # Take 3 month rolling average of monthly rainfall
                     ][,list(Rain.M.Sum3.Mean=mean(Rain.M.Sum3,na.rm=T)),by=list(Index,Month) # Average monthly rainfall by month and site across all years
                       ][,Rain.Season:=SOS_MaxMonth(Rain=Rain.M.Sum3.Mean,Month=Month,Pad=2),by=list(Index) # Code months for 2 x 3 month seasons with highest rainfall
                         ][Rain.Season==2,Rain.Filter:=max(Rain.M.Sum3.Mean)>=Min.Rain,by=list(Index,Rain.Season) # Is season 2  rainfall higher than threshold specified?
                           ][Rain.Filter==F & Rain.Season==2,Rain.Season:=NA] # Set season to NA if rainfall threshold is not met
  
  # Add season code to daily climate information
  CLIM<-merge(CLIM,CLIM.MONTH[,list(Index,Month,Rain.Season)],by=c("Index","Month"),all.x=T)
  
  # Ensure dataset is correctly ordered
  CLIM<-CLIM[order(Index,Year,Dekad)]
  
  # 4) Chunk data for parallel processing #####
  print(paste("Chunking Data:",COUNTRY,"-", match(COUNTRY,Countries),"/",length(Countries)))
  
  # Split list into chunks again based on the number of cores for parallel processing
  ChunkSize<-ceiling(CLIM[,length(unique(Index))]/Cores)
  
  Chunks<-split(CLIM[,unique(Index)], ceiling(seq_along(CLIM[,unique(Index)])/ChunkSize))
  
  Chunks<-rbindlist(lapply(1:length(Chunks),FUN=function(i){
    data.table(Chunk=i,Index=Chunks[[i]])
  }))
  
  CLIM<-merge(CLIM,Chunks,by="Index")
  
  CLIM<-split(CLIM,by="Chunk")
  
  print(paste("Processing:",COUNTRY,"-", match(COUNTRY,Countries),"/",length(Countries)))
  
  # 5) Classify rainy seasons #####
  
  if(Debug==F){
    future::plan(multisession, workers = Cores)
    
    SOS_Data<-future.apply::future_lapply(CLIM,
                                SOS_Wrap,
                                D1.mm=D1.mm,
                                D2.mm=D2.mm,
                                D2.len=D2.len,
                                AI.t=AI.t,
                                PadBack=PadBack,
                                PadForward=PadForward,
                                Do.SeqMerge=Do.SeqMerge,
                                MaxGap=MaxGap,
                                MinStartLen=MinStartLen,
                                MaxStartSep=MaxStartSep,
                                ClipAI=ClipAI,
                                Season2.Prop=Season2.Prop,
                                MinLength=MinLength,
                                AI_Seasonal=AI_Seasonal,
                                RollBack=RollBack,
                                S1.AI=S1.AI,
                                future.packages="terra")
    
    }else{
    
    # Debugging
      SOS_Data<-lapply(1:length(CLIM),FUN=function(i){
        print(paste0("Country: ",COUNTRY," - ",i,"/",length(CLIM)))
        X<-CLIM[[i]]
        SOS_Wrap(DATA=X,  
                 D1.mm=D1.mm,
                 D2.mm=D2.mm,
                 D2.len=D2.len,
                 AI.t=AI.t,
                 PadBack=PadBack,
                 PadForward=PadForward,
                 Do.SeqMerge=Do.SeqMerge,
                 MaxGap=MaxGap,
                 MinStartLen=MinStartLen,
                 MaxStartSep=MaxStartSep,
                 ClipAI=ClipAI,
                 Season2.Prop=Season2.Prop,
                 MinLength=MinLength,
                 AI_Seasonal=AI_Seasonal,
                 RollBack=RollBack,
                 S1.AI=S1.AI)
      }) 
    
      }
    
  # 6) Combine results into a list and save #####
    SOS_Data<-list(
      Dekadal_SOS=rbindlist(lapply(SOS_Data,"[[","Dekadal_SOS")),
      Seasonal_SOS2=rbindlist(lapply(SOS_Data,"[[","Seasonal_SOS2")),
      LTAvg_SOS2=rbindlist(lapply(SOS_Data,"[[","LTAvg_SOS2")),
      Seasonal_SOS3=rbindlist(lapply(SOS_Data,"[[","Seasonal_SOS3")),
      LTAvg_SOS3=rbindlist(lapply(SOS_Data,"[[","LTAvg_SOS3")),
      LTAVg_Data=CLIM.LT,
      LTAVg_Summary=CLIM.LT.Summary
      )
    
  
    print(paste("Saving:",COUNTRY,"-", match(COUNTRY,Countries),"/",length(Countries)))
    
    save(SOS_Data,file=FILE,compress="gzip",compression_level=6)
    
    rm(SOS_Data,CLIM,CLIM.LT,CLIM.LT.Summary)
    gc()
  }
}
