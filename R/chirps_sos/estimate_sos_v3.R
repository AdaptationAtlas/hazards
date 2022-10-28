# https://github.com/tidyverse/lubridate/issues/617#issuecomment-521393447
#' Format date as yearly or monthly dekad (10-day period)

require(lubridate)
require(data.table)
require(ggplot2)
require(miceadds)
require(sp)
require(pbapply)
require(doSNOW)
require(terra)
require(raster)
require(circular)


source("chirps_sos/sos_functions.R")

options(scipen=999)

# 1) Read in and processes a time-series of daily climate data ####

# Dataset requires fields: 
# 1) Index (Unique identity for location); 
# 2) Date;
# 3) Rain (daily rainfall in mm);
# 4) ETo (potential evapotranspiration in mm).

# Data and save locations  ####

CHIRPS_Dir<-"C:/Datasets/CHIRPS/CountryLists"
Hobbins_Dir<-"C:/Datasets/HobbinsRefET/Countries"
CHIRPS_Dekad_Dir<-paste0("C:/Datasets/CHIRPS/CountryLists","/Dekad")
CHIRPS_SOS_Dir<-"C:/Datasets/CHIRPS/SOS/Countries"

if(!dir.exists(CHIRPS_SOS_Dir)){
  dir.create(CHIRPS_SOS_Dir,recursive = T)
}

if(!dir.exists(CHIRPS_Dekad_Dir)){
  dir.create(CHIRPS_Dekad_Dir)
}

CHIRPS_Ref<-load.Rdata2(list.files(CHIRPS_Dir,"CellReference"),path=CHIRPS_Dir)
CHIRPS_Times2<-load.Rdata2(list.files(CHIRPS_Dir,"Times"),path=CHIRPS_Dir)[,Dekad:=SOS_Dekad(Date,type="year")
                                                                          ][,Year:=as.numeric(format(Date,"%Y"))
                                                                            ][,Month:=as.numeric(format(Date,"%m"))]
CHIRPS_Times<-unique(CHIRPS_Times2[,list(Dekad,Month,Year)])

#Hobbins_Ref<-load.Rdata2(list.files(Hobbins_Dir,"CellReference"),path=Hobbins_Dir)
Hobbins_Times<-load.Rdata2(list.files(Hobbins_Dir,"Times"),path=Hobbins_Dir)

# Load map of Africa
AfricaMap<-rworldmap::getMap(resolution = "high")
AfricaMap<-AfricaMap[AfricaMap$REGION=="Africa"&!is.na(AfricaMap$REGION),]
AfricaMap<-spTransform(AfricaMap,CRSobj="+init=epsg:4326")


# Create Functions ####

# This function takes two angles and calculates the difference between then as such we need to convert dekads 1:36 to degrees 1:360 (and back again)
CicularDist<-function(Val1,Val2){
  as.numeric(circular::dist.circular(circular::circular(c(Val1,Val2),type="angles",units = "degrees"),method="geodesic")*180/pi)
}

SOS_Fun<-function(DATA,
                  D1.mm=25,
                  D2.mm=20,
                  D2.len=2,
                  AI.t=0.5,
                  Do.SeqMerge=T,
                  PadBack=3,
                  PadForward=3,
                  MaxGap=1,
                  MinStartLen=2,
                  MaxStartSep=1,
                  ClipAI=F,
                  AI_Seasonal=F,
                  Skip2=F,
                  S1.AI=T
                  ){
  if(Skip2==F){
  DATA<-DATA[,list(Rain.Dekad=sum(Rain),AI=mean(AI)),by=list(Index,Year,Dekad,Rain.Season) # Sum rainfall and take mean aridity, Index by dekad  (within year, site and rain season)
  ][,Dekad.Season:=SOS_SeasonPad(Data=Rain.Season,PadBack=PadBack,PadForward=PadForward),by=Index] # Pad rainy seasons (for growing season > 150 days)
  }
  
  DATA<-DATA[,Dekad.Seq:=SOS_UniqueSeq(Dekad.Season),by=Index # Sequences within sites need a unique ID
  ][,Complete:=length(Dekad)==36,by=list(Index,Year) # Calculate dekads within a year
  ][Complete==T | (Dekad %in% 34:36 & Year==min(Year)) # Remove incomplete years but keep last three dekads (when wet period start is Jan we need to look 3 dekads before this) 
  ][,Complete:=NULL # Tidy up
  ][,Rain.sum2:=slide_apply(Rain.Dekad,window=D2.len+1,step=1,fun=sum) # Rainfall for next two dekads
  ][,SOSmet:=Rain.sum2>=D2.mm & Rain.Dekad>=D1.mm] # Is rainfall of current dekad >=25 and sum of next 2 dekads >=20?
    
  if(AI_Seasonal==T){
    DATA<-DATA[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad,Year)] # Calculate mean aridity Site.Key per dekad across timeseries
  }else{
    DATA<-DATA[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad)] # Calculate mean aridity Site.Key per dekad across timeseries
  }
  
  DATA<-DATA[,AI.0.5:=AI.mean>=AI.t,by=AI.mean # Is aridity >=0.5?
  ][!(is.na(Dekad.Season)),AI.Seq1:=SOS_RSeason(RAIN=SOSmet,AI=AI.0.5,S1.AI=S1.AI),by=list(Index,Dekad.Seq)] # Look for sequences of AI>=0.5 starting when rainfall criteria met
  
  if(Do.SeqMerge){
    DATA[!(is.na(Dekad.Season)),AI.Seq:=SOS_SeqMerge(Seq=AI.Seq1,AI=AI.0.5,MaxGap=MaxGap,MinStartLen=MinStartLen,MaxStartSep=MaxStartSep,ClipAI=ClipAI,S1.AI=S1.AI),by=list(Index,Dekad.Seq)]
  }else{
    DATA[,AI.Seq:=AI.Seq1]
  }
  
  DATA<-DATA[!is.na(AI.Seq),SOS:=Dekad[1],by=list(Index,AI.Seq,Dekad.Seq) # Start of season (SOS) is first dekad of each sequence
  ][!is.na(AI.Seq),EOS:=Dekad[length(Dekad)],by=list(Index,AI.Seq,Dekad.Seq) # End of season (EOS) is last dekad of each sequence
  ][SOS<EOS,LGP:=EOS-SOS # Length of growing period (LGP) is SOS less EOS
  ][SOS>EOS,LGP:=36-SOS+EOS # Deal with scenario where SOS is in different year to EOS
  ][SOS==EOS,c("AI.Seq","SOS","EOS"):=NA # Remove observations where SOS == EOS (sequence is length 0)
  ][Year==max(Year) & EOS==36,c("LGP","EOS"):=NA # remove EOS and LGP where EOS is the last dekad of the available data
  ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Start.Year:=Year[1],by=list(Index,Dekad.Seq) # Add starting year for seasons
  ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq)] # Add total rainfall for season
  
  return(DATA)
  
}

# Create a wrapper for data.table operations  ####

SOS_Wrap<-function(DATA,
                   D1.mm=25,
                   D2.mm=20,
                   D2.len=2,
                   AI.t=0.5,
                   PadBack=PadBack,
                   PadForward=PadForward,
                   Do.SeqMerge=T,
                   MaxGap=1,
                   MinStartLen=2,
                   MaxStartSep=1,
                   ClipAI=F,
                   Season2.Prop=0.33,
                   MinLength=4,
                   AI_Seasonal=F,
                   RollBack=F,
                   S1.AI=T){
  
  # 1) First pass analysis ####
  
  CLIM.Dekad<-SOS_Fun(DATA,
                      D1.mm=D1.mm,
                      D2.mm=D2.mm,
                      D2.len=D2.len,
                      AI.t=AI.t,
                      Do.SeqMerge=Do.SeqMerge,
                      PadBack=PadBack,
                      PadForward=PadForward,
                      MaxGap=MaxGap,
                      MinStartLen=MinStartLen,
                      MaxStartSep=MaxStartSep,
                      ClipAI=ClipAI,
                      AI_Seasonal = AI_Seasonal,
                      Skip2 = F,
                      S1.AI=S1.AI)
  
  # 2) Calculate Seasonal Values ####
  Len<-CLIM.Dekad[,length(unique(Year))]
  Seasonal<-unique(CLIM.Dekad[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,SOS,EOS,LGP,Dekad.Season,Tot.Rain)])
  Seasonal[!is.na(Dekad.Season),Seasons.Count:=.N,by=list(Index,Dekad.Season)
  ][,Season2Prop:=Seasons.Count/Len,by=Index]
  

  Seasonal[,Seasons:=length(unique(Dekad.Season)),by=Index]
  
  # 3) Roll back SOS where SOS is fixed ####
  if(RollBack==T){
  # Add similarity field (proportion of SOS dekads which are the same as the most frequent SOS dekads)
  SameSOS<-function(SOS){
    N<-length(SOS)
    SOS<-SOS[!is.na(SOS)]
    if(length(SOS)>0){
      X<-table(SOS)
      return(round(max(X)/N,2))
    }else{
      return(NA)
    }
  }
  
  Seasonal[,SOSsimilarity:=SameSOS(SOS),by=list(Index,Dekad.Season)
           ][,SOSNA:=sum(is.na(SOS))/.N,by=list(Index,Dekad.Season)]
  
  
  # Subset to very similar planting dates and sites where NAs are not frequent
  X<-unique(Seasonal[SOSsimilarity>0.95 & SOSNA<0.2,list(Index,Dekad.Season,Seasons)])
  }
  # 3.1) Scenario 1: SOS fixed and one season present #####
  # This is a simple case of rolling back the one season
  if(RollBack==T){
  Sites<-X[Seasons==1,Index]
  
  if(length(Sites)>0){
    # Double padding rainy of season start date
    CLIM.Dekad2<-SOS_Fun(DATA[Index %in% Sites],D1.mm=D1.mm,
                        D2.mm=D2.mm,
                        D2.len=D2.len,
                        AI.t=AI.t,
                        Do.SeqMerge=Do.SeqMerge,
                        PadBack=PadBack*2,
                        PadForward=0,
                        MaxGap=MaxGap,
                        MinStartLen=MinStartLen,
                        MaxStartSep=MaxStartSep,
                        ClipAI=ClipAI,
                        AI_Seasonal = AI_Seasonal,
                        Skip2=F,
                        S1.AI=S1.AI)
    
    CLIM.Dekad<-rbind(CLIM.Dekad[!Index %in% Sites],CLIM.Dekad2)
  }
  }
  # 3.2) Scenario 2: SOS fixed and two seasons present #####
  if(RollBack==T){
  # Subset data
  Data<-CLIM.Dekad[Index %in% X[Seasons==2,Index]]
  
  # Calculate season separation
  SeasonSpacing<-function(SOS,EOS,Dekad.Season){
    if(length(unique(Dekad.Season))>=2){
      Data<-unique(data.table(SOS=SOS,EOS=EOS,Dekad.Season=Dekad.Season))
      
      SOSEOS<-Data[!is.na(Dekad.Season),list(SOS=as.numeric(median(SOS,na.rm = T)),EOS=as.numeric(median(EOS,na.rm = T))),by=list(Dekad.Season)]
      
      # Difference between start season two and end season one 
      SOS<-SOSEOS[Dekad.Season==2,SOS]
      EOS<-SOSEOS[Dekad.Season==1,EOS]
      if(SOS<EOS){
        SOS<-36-EOS+1
        EOS<-1
      }
      Diff.1vs2<-SOS-EOS
      
      # Difference between start season one and end season two
      SOS<-SOSEOS[Dekad.Season==1,SOS]
      EOS<-SOSEOS[Dekad.Season==2,EOS]
      if(SOS<EOS){
        SOS<-36-EOS+1
        EOS<-1
      }
      
      Diff.2vs1<-SOS-EOS
      
      Diffs<-c(Diff.1vs2,Diff.2vs1)
      
      Diff<-data.table(sepmin=min(Diffs)[1],sepmax=max(Diffs)[1],order=which(Diffs==min(Diffs))[1])
      
      
    }else{
      Diff<-data.table(sepmin=as.numeric(NA),sepmax=as.numeric(NA),order=as.numeric(NA))
    }
    
    return(Diff)
  }
  
  # If order ==1 then adjacent seasons are ordered 1 then 2, if 2 vice versa
  Data[,Season.Sep.Min:=SeasonSpacing(SOS,EOS,Dekad.Season)$sepmin,by=Index
  ][,Season.Sep.Max:=SeasonSpacing(SOS,EOS,Dekad.Season)$sepmax,by=Index
  ][,Season.Order:=SeasonSpacing(SOS,EOS,Dekad.Season)$order,by=Index]
  
  # Calculate number of fixed season and flag which seasons are fixed
  X.Seasons<-X[Seasons==2,list(FixedSeasons.N=length(unique(Dekad.Season)),FixedSeasons=paste(unique(Dekad.Season),collapse = "-")),by=Index]
  Data<-merge(Data,X.Seasons,by="Index")
  
  # 3.2.1) Seasons are Adjacent ######
  DataAdjacent<-Data[Season.Sep.Min<2]
  
  # If leading season!=fixed season then there is nothing to change (fixed season immediately adjacent to leading season)
  DataAdjacentFixed1<-DataAdjacent[(FixedSeasons.N==1 & Season.Order==FixedSeasons)|FixedSeasons.N==2]
  
  # If we have adjacent seasons and the first season is fixed (i.e. date need adjusting back) then we adjust the start of the season
  # window for both season 1 and season 2. This should help balance the lengths of the two seasons where the rainy season is long enough
  # to accommodate two growing seasons.
  
  Sites<-DataAdjacentFixed1[,unique(Index)]
  
  if(length(Sites)>0){
    # Double padding of rainy season start date remove padding of end date
    # Note for non-adjacent seasons a different method is used that can accommodate flexible padding length by Site
    CLIM.Dekad1<-SOS_Fun(DATA[Index %in% Sites],D1.mm=D1.mm,
                        D2.mm=D2.mm,
                        D2.len=D2.len,
                        AI.t=AI.t,
                        Do.SeqMerge=Do.SeqMerge,
                        PadBack=PadBack*2,
                        PadForward=0,
                        MaxGap=MaxGap,
                        MinStartLen=MinStartLen,
                        MaxStartSep=MaxStartSep,
                        ClipAI=ClipAI,
                        AI_Seasonal = AI_Seasonal,
                        Skip2=F,
                        S1.AI=S1.AI
                        )
    
    CLIM.Dekad<-rbind(CLIM.Dekad[!Index %in% c(Sites)],CLIM.Dekad1)
    Clim.Dekad1<-NULL
  }
  
  }
  # 3.2.1) Seasons are not adjacent ######
  # Need to count back for fixed season but not beyond EOS of other season
  if(RollBack==T){
  # Subset to seasons with a separation of at least 1 dekad 
  DataNonAdjacent<-Data[Season.Sep.Min>=2] # >= 2 is correct
  
  DataNonAdjacent[Season.Order==2 & Rain.Season==2 & grepl(2,FixedSeasons),Season.Sep:=Season.Sep.Max]
  DataNonAdjacent[Season.Order==2 & Rain.Season==1 & grepl(1,FixedSeasons),Season.Sep:=Season.Sep.Min]
  DataNonAdjacent[Season.Order==1 & Rain.Season==1 & grepl(1,FixedSeasons),Season.Sep:=Season.Sep.Min]
  DataNonAdjacent[Season.Order==1 & Rain.Season==2 & grepl(2,FixedSeasons),Season.Sep:=Season.Sep.Max]
  DataNonAdjacent[!is.na(Rain.Season) & is.na(Season.Sep),Season.Sep:=0]
  
  # Set a limit on maximum number of dekads to roll back           
  DataNonAdjacent[Season.Sep>PadBack,Season.Sep:=PadBack]
  
  Sites<-DataNonAdjacent[,unique(Index)]
  
  if(length(Sites)>0){
    CLIM.Dekad1<-DATA[Index %in% Sites,list(Rain.Dekad=sum(Rain),AI=mean(AI)),by=list(Index,Year,Dekad,Rain.Season)]
    
    # Merge season separation with climate data
    CLIM.Dekad1<-merge(CLIM.Dekad1,
                       unique(DataNonAdjacent[!is.na(Rain.Season),list(Index,Rain.Season,Season.Sep)]),
                       by=c("Index","Rain.Season"),all.x=T)
    
      # Merge fixed season identity with climate data
      CLIM.Dekad1<-merge(CLIM.Dekad1,
                         unique(DataNonAdjacent[!is.na(Rain.Season),list(Index,FixedSeasons)]),
                         by=c("Index"),all.x=T)
    
    # Revert to original order
    CLIM.Dekad1<-CLIM.Dekad1[order(Index,Year,Dekad)]
    
    # Increase padding rainy of season start date and reduce padding of end date
    CLIM.Dekad1[Rain.Season==1,Season1:=Rain.Season
                ][Rain.Season==2,Season2:=Rain.Season
                  ][,Dekad.Season1:=SOS_SeasonPad(Data=Season1,
                                                  PadBack=PadBack+Season.Sep[Rain.Season==1 & !is.na(Rain.Season)][1],
                                                  PadForward=PadForward-Season.Sep[Rain.Season==1 & !is.na(Rain.Season)][1]),by=Index 
                    ][,Dekad.Season2:=SOS_SeasonPad(Data=Season2,
                                                    PadBack=PadBack+Season.Sep[Rain.Season==2 & !is.na(Rain.Season)][1],
                                                    PadForward=PadForward-Season.Sep[Rain.Season==2 & !is.na(Rain.Season)][1]),by=Index 
    ]
    
    # Recombine dekad season numbering
    CLIM.Dekad1[FixedSeasons==1,Dekad.Season:=Dekad.Season1
    ][is.na(Dekad.Season)  &  FixedSeasons==1,Dekad.Season:=Dekad.Season2
    ][FixedSeasons==2,Dekad.Season:=Dekad.Season2
    ][is.na(Dekad.Season)  &  FixedSeasons==2,Dekad.Season:=Dekad.Season1
    ][!FixedSeasons %in% c(1,2),Dekad.Season:=Dekad.Season1
    ][is.na(Dekad.Season)  & !FixedSeasons %in% c(1,2),Dekad.Season:=Dekad.Season2
    ][,Dekad.Season1:=NULL
    ][,Dekad.Season2:=NULL
    ][,Season1:=NULL
    ][,Season2:=NULL
    ][,FixedSeasons:=NULL
    ][,Season.Sep:=NULL]
    
    
    CLIM.Dekad1<-SOS_Fun(DATA=CLIM.Dekad1,
                        D1.mm=D1.mm,
                        D2.mm=D2.mm,
                        D2.len=D2.len,
                        AI.t=AI.t,
                        Do.SeqMerge=Do.SeqMerge,
                        PadBack=PadBack,
                        PadForward=PadForward,
                        MaxGap=MaxGap,
                        MinStartLen=MinStartLen,
                        MaxStartSep=MaxStartSep,
                        ClipAI=ClipAI,
                        AI_Seasonal = AI_Seasonal,
                        Skip2 = T,
                        S1.AI=S1.AI)
    
    if(F){
    CLIM.Dekad1<-CLIM.Dekad1[,Dekad.Seq:=SOS_UniqueSeq(Dekad.Season),by=Index # Sequences within sites need a unique ID
    ][,Complete:=length(Dekad)==36,by=list(Index,Year) # Calculate dekads within a year
    ][Complete==T | (Dekad %in% 34:36 & Year==min(Year)) # Remove incomplete years but keep last three dekads (when wet period start is Jan we need to look 3 dekads before this) 
    ][,Complete:=NULL # Tidy up
    ][,Rain.sum2:=slide_apply(Rain.Dekad,window=D2.len+1,step=1,fun=sum) # Rainfall for next two dekads
    ][,SOSmet:=Rain.sum2>=D2.mm & Rain.Dekad>=D1.mm] # Is rainfall of current dekad >=25 and sum of next 2 dekads >=20?

    if(AI_Seasonal==T){
      CLIM.Dekad1<-CLIM.Dekad1[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad,Year)] # Calculate mean aridity Site.Key per dekad across timeseries
    }else{
      CLIM.Dekad1<-CLIM.Dekad1[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad)] # Calculate mean aridity Site.Key per dekad across timeseries
    }
      
    CLIM.Dekad1<-CLIM.Dekad1[,AI.0.5:=AI.mean>=AI.t,by=AI.mean # Is aridity Index >=0.5?
    ][!(is.na(Dekad.Season)),AI.Seq1:=SOS_RSeason(RAIN=SOSmet,AI=AI.0.5,S1.AI=S1.AI),by=list(Index,Dekad.Seq)] # Look for sequences of AI>=0.5 starting when rainfall criteria met

      if(Do.SeqMerge){
        CLIM.Dekad1[!(is.na(Dekad.Season)),AI.Seq:=SOS_SeqMerge(Seq=AI.Seq1,AI=AI.0.5,MaxGap=MaxGap,MinStartLen=MinStartLen,MaxStartSep=MaxStartSep,ClipAI=ClipAI,S1.AI=S1.AI),by=list(Index,Dekad.Seq)]
      }else{
        CLIM.Dekad1[,AI.Seq:=AI.Seq1]
      }
      
    CLIM.Dekad1<-CLIM.Dekad1[!is.na(AI.Seq),SOS:=Dekad[1],by=list(Index,AI.Seq,Dekad.Seq) # Start of season (SOS) is first dekad of each sequence
    ][!is.na(AI.Seq),EOS:=Dekad[length(Dekad)],by=list(Index,AI.Seq,Dekad.Seq) # End of season (EOS) is last dekad of each sequence
    ][SOS<EOS,LGP:=EOS-SOS # Length of growing period (LGP) is SOS less EOS
    ][SOS>EOS,LGP:=36-SOS+EOS # Deal with scenario where SOS is in different year to EOS
    ][SOS==EOS,c("AI.Seq","SOS","EOS"):=NA # Remove observations where SOS == EOS (sequence is length 1)
    ][Year==max(Year) & EOS==36,c("LGP","EOS"):=NA # remove EOS and LGP where EOS is the last dekad of the available data
    ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Start.Year:=Year[1],by=list(Index,Dekad.Seq) # Add starting year for seasons
    ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq)] # Add total rainfall for season
    }
    
    CLIM.Dekad<-rbind(CLIM.Dekad[!Index %in% Sites],CLIM.Dekad1)
    Clim.Dekad1<-NULL
  }
  }
  # 4.3) Calculate seasonal values #####
  
  Seasonal2<-unique(CLIM.Dekad[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,SOS,EOS,LGP,Dekad.Season,Tot.Rain)])
  # Remove second seasons that are too short
  Seasonal2<-Seasonal2[!(Dekad.Season==2 & LGP<MinLength)]
  Seasonal2<-Seasonal2[!is.na(Dekad.Season),Seasons.Count:=.N,by=list(Index,Dekad.Season)][,Season2Prop:=Seasons.Count/Len]
  
  # Remove second seasons that are present for less than 1/3 the time of first seasons
  if(!is.na(Season2.Prop)){
    Seasonal2<-Seasonal2[Season2Prop>Season2.Prop]
  }
  
  # How many seasons present at a site?
  Seasonal2[,Seasons:=length(unique(Dekad.Season)),by=Index]
  
  # What is the similarity of SOS within the site?
  Seasonal2[,SOSsimilarity:=SameSOS(SOS),by=list(Index,Dekad.Season)]
  
  X1<-unique(Seasonal2[SOSsimilarity>0.95,list(Index,Dekad.Season,Seasons)])

  # 4.4) Add separation #####
  CLIM.Dekad[!(is.na(EOS)|is.na(SOS)|is.na(Dekad.Season)),Season.Sep.Min:=SeasonSpacing(SOS,EOS,Dekad.Season)$sepmin,by=Index
  ][!(is.na(EOS)|is.na(SOS)|is.na(Dekad.Season)),Season.Sep.Max:=SeasonSpacing(SOS,EOS,Dekad.Season)$sepmax,by=Index
  ][!(is.na(EOS)|is.na(SOS)|is.na(Dekad.Season)),Season.Order:=SeasonSpacing(SOS,EOS,Dekad.Season)$order,by=Index]
  
  # 5) Is planting possible in the off season - is this a humid region? ####
  
  # Consider using AIseq here rather than Dekad.Season?
  Sites<-Seasonal2[Seasons==2,unique(Index)]
  
  if(length(Sites)>0){
  CLIM.Dekad1<-data.table::copy(CLIM.Dekad)[Index %in% Sites
  ][is.na(Dekad.Season),Dekad.Season1:=3
  ][!is.na(Dekad.Season),Dekad.Season:=NA
  ][,Dekad.Season:=Dekad.Season1
  ][,Dekad.Season1:=NULL
  ][,Dekad.Seq:=SOS_UniqueSeq(Dekad.Season),by=Index]
  
  # Function to shrink third season by one dekad at each end
  ShrinkX<-function(X){
    X<-unlist(X)
    X[1]<-NA
    X[length(X)]<-NA
    return(X)
  }
  
  if(F){
    CLIM.Dekad1<-CLIM.Dekad1[,Dekad.Seq2:=Dekad.Seq
    ][,Dekad.Seq2:=ShrinkX(Dekad.Seq2),by=list(Index,Dekad.Seq)
    ][,Dekad.Seq:=Dekad.Seq2
    ][,Dekad.Seq2:=NULL
    ]
  }
  
  CLIM.Dekad1<-CLIM.Dekad1[Dekad.Season==3
  ][,Rain.sum2:=slide_apply(Rain.Dekad,window=D2.len+1,step=1,fun=sum) # Rainfall for next two dekads
  ][,SOSmet:=Rain.sum2>=D2.mm & Rain.Dekad>=D1.mm] # Is rainfall of current dekad >=25 and sum of next 2 dekads >=20?
  
    if(AI_Seasonal==T){
      CLIM.Dekad1<-CLIM.Dekad1[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad,Year)] # Calculate mean aridity Site.Key per dekad across timeseries
    }else{
      CLIM.Dekad1<-CLIM.Dekad1[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad)] # Calculate mean aridity Site.Key per dekad across timeseries
    }
  
  CLIM.Dekad1<-CLIM.Dekad1[,AI.0.5:=AI.mean>=AI.t,by=AI.mean # Is aridity Index >=0.5?
  ][!(is.na(Dekad.Season)),AI.Seq1:=SOS_RSeason(RAIN=SOSmet,AI=AI.0.5,S1.AI=S1.AI),by=list(Index,Dekad.Seq)] # Look for sequences of AI>=0.5 starting when rainfall criteria met
  
  if(Do.SeqMerge){
    CLIM.Dekad1[!(is.na(Dekad.Season)),AI.Seq:=SOS_SeqMerge(Seq=AI.Seq1,AI=AI.0.5,MaxGap=MaxGap,MinStartLen=MinStartLen,MaxStartSep=MaxStartSep,ClipAI=ClipAI,S1.AI=S1.AI),by=list(Index,Dekad.Seq)]
  }else{
    CLIM.Dekad1[,AI.Seq:=AI.Seq1]
  }
  
  CLIM.Dekad1<-CLIM.Dekad1[!is.na(AI.Seq),SOS:=Dekad[1],by=list(Index,AI.Seq,Dekad.Seq) # Start of season (SOS) is first dekad of each sequence
  ][!is.na(AI.Seq),EOS:=Dekad[length(Dekad)],by=list(Index,AI.Seq,Dekad.Seq) # End of season (EOS) is last dekad of each sequence
  ][SOS<EOS,LGP:=EOS-SOS # Length of growing period (LGP) is SOS less EOS
  ][SOS>EOS,LGP:=36-SOS+EOS # Deal with scenario where SOS is in different year to EOS
  ][SOS==EOS,c("AI.Seq","SOS","EOS"):=NA # Remove observations where SOS == EOS (sequence is length 1)
  ][Year==max(Year) & EOS==36,c("LGP","EOS"):=NA # remove EOS and LGP where EOS is the last dekad of the available data
  ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Start.Year:=Year[1],by=list(Index,Dekad.Seq) # Add starting year for seasons
  ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq)] # Add total rainfall for season
  
  Seasonal3<-unique(CLIM.Dekad1[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,SOS,EOS,LGP,Dekad.Season,Tot.Rain)])
  # Remove second seasons that are too short
  Seasonal3<-Seasonal3[!(Dekad.Season==3 & LGP<MinLength)]
  Seasonal3[Dekad.Season==3,Seasons.Count:=sum(Dekad.Season==3),by=Index
  ][,Season3Prop:=Seasons.Count/Len,by=Index]
  
  # Remove third seasons that are present for less than 1/3 of the time 
  if(!is.na(Season2.Prop)){
    Seasonal3<-Seasonal3[Season3Prop>Season2.Prop]
  }

  if(nrow(Seasonal3)>0){
    Sites<-Seasonal3[,unique(Index)]
    
    # Combine data main dataset with modified season 3 data
    CLIM.Dekad.3<-rbind(CLIM.Dekad1[Index %in% Sites],CLIM.Dekad[!(Index %in% Sites & is.na(Dekad.Season))])
    
    # Update seasonal statistics
    Seasonal3<-unique(CLIM.Dekad.3[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,Dekad.Season,SOS,EOS,LGP,Tot.Rain)])
    Seasonal3<-Seasonal3[base::order(Index,Start.Year,Dekad.Season,decreasing=c(FALSE,FALSE,FALSE),method="radix")]
    # Remove second seasons that are too short
    Seasonal3<-Seasonal3[!(Dekad.Season %in% c(2,3) & LGP<MinLength)]
    Seasonal3[,Seasons.Count:=.N,by=list(Index,Dekad.Season)
    ][,SeasonProp:=Seasons.Count/Len,by=Index]
    
    # Remove third seasons that are present for less than a specified proportion the time (relative to season 1)
    if(!is.na(Season2.Prop)){
    Seasonal3<-Seasonal3[SeasonProp>=Season2.Prop]
    }
    
    Seasonal3[,Seasons:=length(unique(Dekad.Season)),by=Index]
    
    Seasonal3[,SOSsimilarity:=SameSOS(SOS),by=list(Index,Dekad.Season)]
  }
  }else{
    Seasonal3<-Seasonal2[0]
  }
  
  # 6) Site.Details ####
  # Proportion of dekads, entire time series, where SOS or AI rule is true
  Site.Details<-CLIM.Dekad[,list(AI.0.5.Prop=sum(AI.0.5)/.N,
                                 AI.1.0.Prop=sum(AI.mean>=1)/.N,
                                 SOSmet.Prop=sum(SOSmet,na.rm = T)/.N,
                                 MAP=sum(Rain.Dekad)/length(unique(Year))),by=Index]
  
  
  # 7) Long term average SOS, EOS, LGP and Total Rainfall ####
  
  LTAvg_SOS2<-Seasonal2[!is.na(Dekad.Season),list(SOS.mean=round(CircMean(m=SOS,interval=36,na.rm=T),1),
                                                  SOS.min=suppressWarnings(min(SOS,na.rm=T)),
                                                  SOS.max=suppressWarnings(max(SOS,na.rm=T)),
                                                  SOS.sd=suppressWarnings(sd(SOS,na.rm=T)),
                                                  Total.Seasons=.N,
                                                  EOS.mean=round(CircMean(m=EOS,interval=36,na.rm=T),1),
                                                  LGP.mean=round(mean(LGP,na.rm=T),1),
                                                  LGP.median=round(median(LGP,na.rm=T),1),
                                                  Tot.Rain.mean=round(mean(Tot.Rain,na.rm=T),1),
                                                  Tot.Rain.sd=round(sd(Tot.Rain,na.rm=T),1),
                                                  SOS.EOS.XYearEnd=round(sum(EOS[!is.na(EOS)]<SOS[!is.na(EOS)])/length(EOS[!is.na(EOS)]),2),
                                                  SOS.add15.XYearEnd=round(sum((SOS[!is.na(SOS)]+15)>36,na.rm=T)/length(SOS[!is.na(SOS)]),2)),
                        by=list(Index,Dekad.Season)]
  
  LTAvg_SOS2[!is.na(Dekad.Season),Seasons:=length(unique(Dekad.Season)),by=Index]
  
  # Order seasons by SOS
  LTAvg_SOS2[!(is.na(SOS.mean)|is.na(Seasons)|is.na(Dekad.Season)),Season.Ordered:=(1:length(Dekad.Season))[order(SOS.mean)],by=Index
  ][Seasons==1,Season.Ordered:=NA]
  
  # Merge LT season order and year end data with Seasonal
  Seasonal2<-merge(Seasonal2,LTAvg_SOS2[,list(Index,Dekad.Season,Season.Ordered,SOS.EOS.XYearEnd,SOS.add15.XYearEnd,SOS.min,SOS.max,Total.Seasons)],by=c("Index","Dekad.Season"),all.x=T)
  
  
  LTAvg_SOS3<-Seasonal3[!is.na(Dekad.Season),list(SOS.mean=round(CircMean(m=SOS,interval=36,na.rm=T),1),
                                                  SOS.median=as.numeric(median(SOS,na.rm=T)),
                                                  SOS.min=suppressWarnings(min(SOS,na.rm=T)),
                                                  SOS.max=suppressWarnings(max(SOS,na.rm=T)),
                                                  SOS.sd=suppressWarnings(sd(SOS,na.rm=T)),
                                                  Total.Seasons=.N,
                                                  EOS.mean=round(CircMean(m=EOS,interval=36,na.rm=T),1),
                                                  LGP.mean=round(mean(LGP,na.rm=T),1),
                                                  LGP.median=round(median(LGP,na.rm=T),1),
                                                  Tot.Rain.mean=round(mean(Tot.Rain,na.rm=T),1),
                                                  Tot.Rain.sd=round(sd(Tot.Rain,na.rm=T),1),
                                                  SOS.EOS.XYearEnd=round(sum(EOS[!is.na(EOS)]<SOS[!is.na(EOS)])/length(EOS[!is.na(EOS)]),2),
                                                  SOS.add15.XYearEnd=round(sum((SOS[!is.na(SOS)]+15)>36,na.rm=T)/length(SOS[!is.na(SOS)]),2)),
                        by=list(Index,Dekad.Season)]
  
  LTAvg_SOS3[!is.na(Dekad.Season),Seasons:=length(unique(Dekad.Season)),by=Index]
  
  # Order seasons by SOS
  LTAvg_SOS3[!(is.na(SOS.mean)|is.na(Seasons)|is.na(Dekad.Season)),Season.Ordered:=(1:length(Dekad.Season))[order(SOS.mean)],by=Index
  ][Seasons==1,Season.Ordered:=NA]
  
  # Merge LT season order and year end data with Seasonal
  Seasonal3<-merge(Seasonal3,LTAvg_SOS3[,list(Index,Dekad.Season,Season.Ordered,SOS.EOS.XYearEnd,SOS.add15.XYearEnd,SOS.min,SOS.max,Total.Seasons)],by=c("Index","Dekad.Season"),all.x=T)
  
  # 8) Save and return data ####
  ERA_SOS<-list(Dekadal_SOS=CLIM.Dekad,
                Seasonal_SOS2=if(nrow(Seasonal2)>0){Seasonal2}else{NULL},
                LTAvg_SOS2=if(nrow(LTAvg_SOS2)>0){LTAvg_SOS2}else{NULL},
                Seasonal_SOS3=if(nrow(Seasonal3)>0){Seasonal3}else{NULL},
                LTAvg_SOS3=if(nrow(LTAvg_SOS3)>0){LTAvg_SOS3}else{NULL})
  return(ERA_SOS)
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
Overwrite<-F # Overwrite previous saved analysis?

S1.AI<-T # Set AI to T when for the first value of sequences of RAIN == T

Countries<-as.character(AfricaMap$ADMIN)
NotInHobbins<-c("Cape Verde","Mauritius","Saint Helena","Seychelles")
Issue<-c("Djibouti","Western Sahara")
Countries<-Countries[!Countries %in% c(NotInHobbins,Issue)]
#Countries<-c("Kenya","Ethiopia","South Africa","Malawi","Nigeria","United Republic of Tanzania")

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

SaveDir<-paste0(CHIRPS_SOS_Dir,"/",FolderName)

if(!dir.exists(SaveDir)){
  dir.create(SaveDir,recursive=T)
}

# Loop SOS analysis over chunks ####
for(COUNTRY in Countries){

  FILE<-paste0(SaveDir,"/SOS_",COUNTRY,".RData")

if(!file.exists(FILE) | Overwrite==T){
print(paste("Loading Country",COUNTRY," - ", match(COUNTRY,Countries),"/",length(Countries)))

HOBBINS<-load.Rdata2(file=paste0(COUNTRY,".RData"),path=Hobbins_Dir)

if(nrow(Hobbins_Times)>nrow(CHIRPS_Times)){
  X<-data.table::copy(Hobbins_Times)[,paste(c(Dekad,Year),collapse=""),by=N][,V1]
  Y<-data.table::copy(CHIRPS_Times)[,paste(c(Dekad,Year),collapse=""),by=list(Year,Month,Dekad)][,V1]
  Z<-X %in% Y
  
  HOBBINS<-lapply(HOBBINS,FUN=function(Data){
    Data[Z]
  })
  
  rm(X,Y,Z)
}

# 1) Prepare Data ####

# Sum CHIRPS rain by dekads

CHIRPS_dekad_file<-paste0(CHIRPS_Dekad_Dir,"/",COUNTRY,"_dekads.RData")
if(!file.exists(CHIRPS_dekad_file)){
  print(paste("Restructuring from Day to Dekad -",COUNTRY))
  CHIRPS<-miceadds::load.Rdata2(paste0(COUNTRY,".RData"),path=CHIRPS_Dir)

Indices<-names(CHIRPS)

CHIRPS<-pblapply(CHIRPS,FUN=function(X){
  X[X==-9999]<-NA
  data.table(Rain=X,CHIRPS_Times2[,list(Year,Dekad)])[,list(Rain=sum(Rain)),by=c("Year","Dekad")][,Rain]
}) 

names(CHIRPS)<-Indices

save(CHIRPS,file=CHIRPS_dekad_file)

}else{
  CHIRPS<-miceadds::load.Rdata2(paste0(COUNTRY,"_dekads.RData"),path=CHIRPS_Dekad_Dir)
}

  # Make a gigantic data.table
CLIM<-data.table(
  ETo=round(unlist(HOBBINS),1),
  Rain=round(unlist(CHIRPS),1),
  Dekad=as.integer(CHIRPS_Times[,rep(Dekad,length(CHIRPS))]),
  Year=as.integer(CHIRPS_Times[,rep(Year,length(CHIRPS))]),
  Month=as.integer(CHIRPS_Times[,rep(Month,length(CHIRPS))]),
  Index=as.integer(rep(names(CHIRPS),each=length(CHIRPS[[1]]))) # ChunkIndex
)

# Remove unneeded large files from memory
rm(CHIRPS,HOBBINS)
gc()

if(rmNAs){
  NANIndexs<-CLIM[,list(NAN=any(is.nan(ETo))),by=Index][NAN==T,Index]
  CLIM<-CLIM[!Index %in% NANIndexs]
  CLIM<-CLIM[!is.na(Rain)]
}

# Calculate aridity index (AI)
CLIM<-CLIM[,AI:=round(Rain/ETo,2)] 

format(object.size(CLIM),units="Gb")

# x) Calculate seasonality using long-term averages only

CLIM.LT<-data.table::copy(CLIM)[,Rain.sum2:=slide_apply(Rain,window=D2.len+1,step=1,fun=sum),by=Index
                                ][,Rain.sum9:=slide_apply2(Rain,window=9,step=1,fun=sum),by=Index # Rainfall over 9 dekads (3 months)
                                    ][,list(ETo=mean(ETo,na.rm=T),
                                            Rain=mean(Rain,na.rm=T),
                                            Rain.sum2=mean(Rain.sum2,na.rm=T),
                                            Rain.sum9=mean(Rain.sum9,na.rm=T)),by=list(Index,Dekad)
                                      ][,AI:=round(Rain/ETo,2)
                                       ][,AImet:=AI>=AI.t
                                         ][,SOSmet:=Rain.sum2>=D2.mm & Rain>=D1.mm
                                          ][,Seq2:=SOS_RSeasonLT(RAIN=SOSmet,AI=AImet,S1.AI=S1.AI),by=Index]


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

# 2) Find wet months ####

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

# 3) Chunk data for parallel processing ####

# Split list into chunks again based on the number of cores for parallel processing
Cores<-14
ChunkSize<-ceiling(CLIM[,length(unique(Index))]/Cores)

Chunks<-split(CLIM[,unique(Index)], ceiling(seq_along(CLIM[,unique(Index)])/ChunkSize))

Chunks<-rbindlist(lapply(1:length(Chunks),FUN=function(i){
  data.table(Chunk=i,Index=Chunks[[i]])
}))

CLIM<-merge(CLIM,Chunks,by="Index")

CLIM<-split(CLIM,by="Chunk")

print(paste("Processing:",COUNTRY,"-", match(COUNTRY,Countries),"/",length(Countries)))

# 4) Classify rainy seasons ####

if(Debug==F){
  cl<-makeCluster(Cores)
  clusterEvalQ(cl, list(library(zoo),library(data.table)))
  clusterExport(cl,list("SOS_SeasonPad","SOS_UniqueSeq","SOS_RSeason","SOS_SeqMerge","SOS_Wrap","PadBack","PadForward",
                        "MinLength","D1.mm","D2.mm","D2.len","AI.t","Do.SeqMerge","MaxGap","MinStartLen","MaxStartSep",
                        "ClipAI","Season2.Prop","MinLength","SOS_Fun","AI_Seasonal","RollBack","CircMean",
                        "SOS_RSeasonLT","LTSeqMerge","OrderDekadSeq","slide_apply","slide_apply2","S1.AI"),envir=environment())
  registerDoSNOW(cl)
 
  system.time( 
  SOS_Data<-parLapply(cl,CLIM,fun=function(X){
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
  )
  
  stopCluster(cl)
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
  
  SOS_Data<-list(
    Dekadal_SOS=rbindlist(lapply(SOS_Data,"[[","Dekadal_SOS")),
    Seasonal_SOS2=rbindlist(lapply(SOS_Data,"[[","Seasonal_SOS2")),
    LTAvg_SOS2=rbindlist(lapply(SOS_Data,"[[","LTAvg_SOS2")),
    Seasonal_SOS3=rbindlist(lapply(SOS_Data,"[[","Seasonal_SOS3")),
    LTAvg_SOS3=rbindlist(lapply(SOS_Data,"[[","LTAvg_SOS3")),
    LTAVg_Data=CLIM.LT,
    LTAVg_Summary=CLIM.LT.Summary
    )
  

  if(Debug){
    print(paste("Saving:",COUNTRY,"-", match(COUNTRY,Countries),"/",length(Countries)))
  }
  
  save(SOS_Data,file=FILE)
  
  rm(SOS_Data,CLIM,CLIM.LT,CLIM.LT.Summary)
  gc()
}
}

# Wrangle SOS chunks ####
SOSFiles<-list.files(SaveDir,".RData",full.names = T)
SOSFiles<-data.table(File=SOSFiles,Country=gsub("[.]RData","",tstrsplit(SOSFiles,"_",keep=2)[[1]]))

AnalysisDir<-paste0("SOS/Analysis")
if(!dir.exists(AnalysisDir)){
  dir.create(AnalysisDir)
}

File<-paste0(AnalysisDir,"/",FolderName,".RData")

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
  
  
  CHIRPS_Ref<-load.Rdata2(list.files(CHIRPS_Dir,"CellReference"),path=CHIRPS_Dir)
  

#***************************************#
# Calculate only from LTAVg_Summary ####
MinLength<-1 # Minimum length of season allowed
Min.Rain<-0 # Minimum seasonal rainfall allowed
SeasonGapT<-3 # if the gap between seasons is <=SeasonGapT then it is considered as a split season (i.e. bimodal)

PlotDir<-paste0("SOS/Plots/",FolderName,"/ML",MinLength,"-MR",Min.Rain,"-SG",SeasonGapT)

if(!dir.exists(PlotDir)){
  dir.create(PlotDir,recursive=T)
}

# The code will modify EOS, SOS and LGP according to SeasonGapT and attempt to tidy seasonal organization between S1 and S2

for(COUNTRY in Countries){
  print(COUNTRY)
  
  Map.Subset<-AfricaMap[AfricaMap$ADMIN==COUNTRY,]
  
  SOSData1<-data.table::copy(SOSlist[[COUNTRY]]$LTAVg_Summary)
  
  SOSData1[,SOS:=round(SOS,0)][SOS==0,SOS:=36]
  SOSData1[,EOS:=round(EOS,0)][EOS==0,EOS:=36]
 
  # NOTE - here seasons with insufficient length or rainfall are removed. This could be done after merging of bimodal seasons?
  SOSData1<-SOSData1[LGP>=MinLength & Tot.Rain>Min.Rain]
  
  suppressWarnings(SOSData1[,Seq:=stringi::stri_replace_all_regex(Seq,
                         pattern=unique(Seq[!is.na(Seq)]),
                         replacement=1:length(unique(Seq[!is.na(Seq)]))),by=Index])
  
  SOSData1<-merge(CHIRPS_Ref,SOSData1,by="Index",all.x=F)

  SOSData1<-dcast(SOSData1,Index+Row+Col+Seasons~Seq,value.var = c("SOS","EOS","LGP","Tot.Rain","Tot.ETo","Balance"))
  if(!"SOS_2" %in% colnames(SOSData1)){
    SOSData1[,c("SOS_2","EOS_2","LGP_2","Tot.Rain_2","Tot.ETo_2","Balance_2"):=NA]
  }
  
  SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist12:=round(CicularDist(EOS_1*10,SOS_2*10)/10,0),by=Index]
  SOSData1[!(is.na(SOS_1)|is.na(SOS_2)),Dist21:=round(CicularDist(EOS_2*10,SOS_1*10)/10,0),by=Index]

  SOSData1[,MinDist:=paste(which(c(Dist12,Dist21)<=SeasonGapT),collapse="-"),by=Index][MinDist=="",MinDist:=NA]
  SOSData1[,MinDistN:=sum(c(Dist12,Dist21)<=SeasonGapT),by=Index]
  
  SOSData<-data.table::copy(SOSData1)
  
  # Consider using S1, S2, etc. to record modified values
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
  
  SOSData[,Seasons:=sum(c(!is.na(S1),!is.na(S2)))*10,by=Index]
  # Manually tidy up countries which are largely unimodal but have a few pixels that split the seasons.
  
  # Check that correct LGP field is being used here
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
          ][Seasons==20,System:="2 WS"
            ][Seasons==10 & MinDistN %in% c(1,2),System:="1 WS-Split"
                 ][Seasons==10 & !MinDistN %in% c(1,2),System:="1 WS"
                  ][LGPsum>=30,System:=paste0(System,"-LGP>=30")
                    ][LGPsum<=6,System:=paste0(System,"-LGP<=6")]
  
  if(!dir.exists(paste0(PlotDir,"/Data"))){
    dir.create(paste0(PlotDir,"/Data"))
  }
  Colnames<-colnames(SOSData)[!grepl("_3|_4|_5|_6",colnames(SOSData))]
  SOSData<-SOSData[,..Colnames]
  save(SOSData,file=paste0(PlotDir,"/Data/",COUNTRY,".RData"))
    
  # prepare coordinates, data, and proj4string
  coords <- data.frame(SOSData[ , c("Col", "Row")])   # coordinates
  crs    <- CRS("+init=epsg:4326") # proj4string of coords
  
  # make the SpatialPointsDataFrame object
  X <- terra::vect(SpatialPointsDataFrame(coords      = coords,
                                          data        = SOSData , 
                                          proj4string = crs))
  
  #BaseRaster<-terra::rast("C:/Datasets/CHIRPS/chirps-v2.0.1981.01.10.tif")
  
  BaseRaster<-terra::rast(ncol=SOSData[,length(unique(Col))], 
                          nrow=SOSData[,length(unique(Row))], 
                          xmin=SOSData[,min(Col)], 
                          xmax=SOSData[,max(Col)], 
                          ymin=SOSData[,min(Row)],
                          ymax=SOSData[,max(Row)])
  
  # SOS Maps ####
  Y1<-round(terra::rasterize(X,BaseRaster,field="S1"),0)
  if(SOSData[!is.na(S2),.N==0]){
    Y2<-Y1
    Y2[]<-NA
  }else{
    Y2<-round(terra::rasterize(X,BaseRaster,field="S2"),0)
  }
  Ymax<-round(terra::rasterize(X,BaseRaster,field="Swettest"),0)
  Yseasons<-terra::rasterize(X,BaseRaster,field="Seasons")
  
  LGP1<-round(terra::rasterize(X,BaseRaster,field="LGP1"),0)
  
  if(!20 %in% SOSData$Seasons){
    LGP2<-LGP1
    LGP2[]<-NA
  }else{
    LGP2<-round(terra::rasterize(X,BaseRaster,field="LGP2"),0)
  }
  
  
  Ystack<-c(Y1,Y2,Ymax,LGP1,LGP2,Yseasons)

  names(Ystack)<-c("SOS_1",
                   "SOS_2",
                   "SOS_Wettest",
                   "LGP_1",
                   "LGP_2",
                   "Seasons")
                   
  
  Ystack2<-raster::stack(Ystack)
  
  coords <- xyFromCell(Ystack2, seq_len(ncell(Ystack2)))
  Data <- raster::stack(as.data.frame(raster::getValues(Ystack2)))
  names(Ystack2) <- names(Ystack)
  
  Data <- cbind(coords, Data)
  Data$values<-factor(Data$values,levels=1:36)
  
  SOSmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=3) +
    #viridis::scale_fill_viridis(option="turbo",discrete = T,drop=F) +
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
  
  # System Maps ####
  
  Ysystem<-terra::rasterize(X,BaseRaster,field="System")
  names(Ysystem)<-"System"
  coords <- xyFromCell(Ysystem, seq_len(ncell(Ysystem)))
  Data <-values(Ysystem)
  names(Ysystem) <- names(Ysystem)
  
  Levels<-data.table(value=sort(unique(Ysystem[!is.na(Ysystem)])),description=levels(Ysystem)[[1]])
  
  
  Data <- data.frame(cbind(coords, Data))
  Data$System<-Levels$description[match(Data$System,Levels$value)]
  
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
  
  g<-gridExtra::arrangeGrob(grobs=list(SOSmap,Ysystem),layout_matrix=matrix(c(rep(1,3),2),nrow = 1))
  
  ggsave(filename = paste0(COUNTRY,".png"),
         plot = g,
         path = PlotDir,
         width= 180,
         height = 100,
         units = "mm",
         scale = 1.2,
         dpi = 600,
         type = "cairo",
         bg="white")
  
}

  # Plot all countries ####
  Files<-list.files(paste0(PlotDir,"/Data"),".RData",full.names = T)
  
  SOSData<-data.table::rbindlist(lapply(Files,miceadds::load.Rdata2),use.names = T)
  SOSData[,Seasons:=sum(!is.na(SOS_1),!is.na(SOS_2))]
  
  # prepare coordinates, data, and proj4string
  coords <- data.frame(SOSData[ , c("Col", "Row")])   # coordinates
  crs    <- CRS("+init=epsg:4326") # proj4string of coords
  
  # make the SpatialPointsDataFrame object
  X <- terra::vect(SpatialPointsDataFrame(coords      = coords,
                                          data        = SOSData , 
                                          proj4string = crs))
  
  BaseRaster<-terra::rast(ncol=SOSData[,length(unique(Col))], 
                          nrow=SOSData[,length(unique(Row))], 
                          xmin=SOSData[,min(Col)], 
                          xmax=SOSData[,max(Col)], 
                          ymin=SOSData[,min(Row)],
                          ymax=SOSData[,max(Row)])
  
  # SOS Maps ####
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
  Ystack<-c(Y1,Y2,Ymax,LGP1,LGP2,Yseasons)
  
  names(Ystack)<-c("SOS_1",
                   "SOS_2",
                   "SOS_Wettest",
                   "LGP_1",
                   "LGP_2",
                   "Seasons")
  
  Ystack2<-raster::stack(Ystack)
  
  coords <- xyFromCell(Ystack2, seq_len(ncell(Ystack2)))
  Data <- raster::stack(as.data.frame(raster::getValues(Ystack2)))
  names(Ystack2) <- names(Ystack)
  
  plot(Ystack)
  
  Data <- cbind(coords, Data)
  Data$values<-factor(Data$values,levels=1:36)
  
  SOSmap<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = values)) +
    facet_wrap(~ ind,ncol=3) +
    #viridis::scale_fill_viridis(option="turbo",discrete = T,drop=F) +
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
         width= 180,
         height = 100,
         units = "mm",
         scale = 1.2,
         dpi = 600,
         type = "cairo",
         bg="white")
  
  # System Map ####
  Ysystem<-terra::rasterize(X,BaseRaster,field="System")
  names(Ysystem)<-"System"
  
  terra::writeRaster(Ysystem,paste0(PlotDir,"/Data/Africa_Systems.tif"),overwrite=T)
  
  
  coords <- xyFromCell(Ysystem, seq_len(ncell(Ysystem)))
  Data <-values(Ysystem)
  names(Ysystem) <- names(Ysystem)
  
  Levels<-data.table(value=sort(unique(Ysystem[!is.na(Ysystem)])),description=levels(Ysystem)[[1]])
  
  
  Data <- data.frame(cbind(coords, Data))
  Data$System<-Levels$description[match(Data$System,Levels$value)]
  
  Ysystem<-ggplot(Data) + 
    geom_tile(aes(x, y, fill = System)) +
    viridis::scale_fill_viridis(option="turbo",discrete = T,na.value="transparent",na.translate = F)+
    theme_bw()+
    labs(fill="System")+
    theme(legend.position = "right",
          axis.title = element_blank(),
          legend.direction = "vertical")+
    guides(fill=guide_legend(ncol=1))+
    coord_equal()
  
  plot(Ysystem)
  
  
  ggsave(filename = "Africa - systems.png",
         plot = Ysystem,
         path = PlotDir,
         width= 180,
         height = 100,
         units = "mm",
         scale = 1.2,
         dpi = 600,
         type = "cairo",
         bg="white")
  
  
  
  
  #*****************************************************#
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
    
    # prepare coordinates, data, and proj4string
    coords <- data.frame(SOSData.S1[ , c("Col", "Row")])   # coordinates
    crs    <- CRS("+init=epsg:4326") # proj4string of coords
    
    # make the SpatialPointsDataFrame object
    X <- terra::vect(SpatialPointsDataFrame(coords      = coords,
                                            data        = SOSData.S1 , 
                                            proj4string = crs))
    
    #BaseRaster<-terra::rast("C:/Datasets/CHIRPS/chirps-v2.0.1981.01.10.tif")
    
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
