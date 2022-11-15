require(lubridate)
require(data.table)
require(miceadds)
require(future)
require(terra)
require(circular)
require(sp)

options(scipen=999)

#Cores<-parallel::detectCores()
Cores<-5

# Data and save locations  ####
DataDir<-"/home/jovyan/common_data"

# Version
version<-1

# Set directories
CHIRPSraw_Dir<-paste0(DataDir,"/chirps_af/raw")
CHIRPS_Dir<-paste0(DataDir,"/chirps_af/intermediate/array_countries")
Hobbins_Dir<-paste0(DataDir,"/hobbins_ref_et/intermediate/array_countries")
CHIRPS_Dekad_Dir<-paste0(DataDir,"/chirps_af/intermediate/array_countries_dekad")
SOS_Dir<-paste0(DataDir,"/atlas_SOS/intermediate/v",version)

if(!dir.exists(SOS_Dir)){
  dir.create(SOS_Dir,recursive = T)
}

if(!dir.exists(CHIRPS_Dekad_Dir)){
  dir.create(CHIRPS_Dekad_Dir)
}

# Import Functions ####
#source("R/chirps_sos/sos_functions.R")
source("sos_functions.R")

# Create a wrapper for data.table operations  ####
#' @param DATA data.frame, data.table or tibble of dekadal climate data. Must contain the fields `Index`, `Dekad`, `Year`, `Rain.Season`,`Rain`
#' `ETo`.
#' @param D1.mm amount of rainfall (mm) that needs to fall in a dekad to trigger onset of rain (dekad_n). Default = 25mm.
#' @param D2.mm amount of rainfall (mm) that needs to fall in subsequent dekads `(dekad_n + 1):(dekad_n + D2.len)`, as defined in `D2.len`, to trigger onset of rain (SOS). Default = 20mm.
#' @param D2.len the number of dekads after dekad 1 from which `D2.mm` is summed. Default = 2.
#' @param AI.t the aridity index (AI) threshold, if AI falls below this thresholds this defines the end of season (EOS). Default = 0.5.
#' @param Do.SeqMerge  logical, if `T` then sequences that are close together are merged, and remove false starts. Default = 1.
#' @param MaxGap an integer value describing the maximum gap (number of NA values) allowed between non-NA values before the sequence breaks.
#' @param MinStartLen an integer value describing the minimum length of first sequence block, if the first sequence is too short and too separated from the next sequence it is removed. This is to remove false starts. Default = 2.
#' @param MaxStartSep an integer value describing the maximum separation of the first sequence block, if the first sequence is too short and too separated from the next sequence it is removed. This is to remove false starts. Default = 1.
#' @param ClipAI logical `T/F`, if `T` then `AI` values corresponding to the last non-NA value in `Seq` are all set to `F` and the sequence is halted by the last `F` value of AI. Default = F.
#' @param S1.AI  a logical vector, when S1.AI==T then the corresponding AI value to the first value in a sequence of RAIN==T is set to `TRUE`. If `S1.AI` is set to `TRUE` in the functions that generate the `Seq` parameter it should also be set to `TRUE` here. Default = T.
#' @param PadBack an integer value specifying number of places to increase sequence backwards. This should not be greater than the minimum separation of sequences. Default = 3.
#' @param PadForward an integer value specifying number of places to increase sequence forwards. This should not be greater than the minimum separation of sequences. Default = 3.
#' @param AI_Seasonal logical, if `T` the seasonal AI is used to define growing seasons, if `F` the long-term average is used. Default = `F`.
#' @param Skip2 logical, if `T`

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
  
  DATA<-data.table(DATA)
  
  if(Skip2==F){
  DATA<-DATA[,list(Rain.Dekad=sum(Rain),AI=mean(AI)),by=list(Index,Year,Dekad,Rain.Season) # Sum rainfall and take mean aridity, Index by dekad  (within year and rain season)
  ][,Dekad.Season:=SOS_SeasonPad(Data=Rain.Season,PadBack=PadBack,PadForward=PadForward),by=Index] # Pad rainy seasons (for growing season > 150 days)
  }
  
  DATA<-DATA[,Dekad.Seq:=SOS_UniqueSeq(Dekad.Season),by=Index # Sequences within sites need a unique ID
  ][,Complete:=length(Dekad)==36,by=list(Index,Year) # Calculate dekads within a year
  ][Complete==T | (Dekad %in% 34:36 & Year==min(Year)) # Remove incomplete years but keep last three dekads (when wet period start is Jan we need to look 3 dekads before this) 
  ][,Complete:=NULL # Tidy up
  ][,Rain.sum2:=slide_apply(Rain.Dekad,window=D2.len+1,step=1,fun=sum) # Rainfall for next two dekads
  ][,SOSmet:=Rain.sum2>=D2.mm & Rain.Dekad>=D1.mm] # Is rainfall of current dekad >=25 and sum of next 2 dekads >=20?
    
  if(AI_Seasonal){
    DATA<-DATA[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad,Year)] # It may be sufficient to simply set AI.mean to AI
  }else{
    DATA<-DATA[,AI.mean:=round(mean(AI,na.rm=T),2),by=list(Index,Dekad)] # Calculate mean aridity index for per dekad across all years
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
#' @param Season2.Prop
#' @param MinLength
#' @param RollBack
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
  
  X1<-unique(Seasonal2[SOSsimilarity>0.95,list(Index,Dekad.Season,Seasons)]) # Redundant?

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
  
  LTAvg_SOS2<-Seasonal2[!is.na(Dekad.Season),list(Total.Seasons=.N,
                                                  SOS.mean=round(CircMean(m=SOS,interval=36,na.rm=T),1),
                                                  SOS.median=as.numeric(median(SOS,na.rm=T)),
                                                  SOS.mode=getmode(SOS,na.rm=T),
                                                  SOS.min=suppressWarnings(min(SOS,na.rm=T)),
                                                  SOS.max=suppressWarnings(max(SOS,na.rm=T)),
                                                  SOS.sd=suppressWarnings(sd(SOS,na.rm=T)),
                                                  EOS.mean=round(CircMean(m=EOS,interval=36,na.rm=T),1),
                                                  EOS.mode=getmode(EOS,na.rm=T),
                                                  EOS.median=round(median(EOS,na.rm=T),1),
                                                  EOS.min=suppressWarnings(min(EOS,na.rm=T)),
                                                  EOS.max=suppressWarnings(max(EOS,na.rm=T)),
                                                  EOS.sd=suppressWarnings(sd(EOS,na.rm=T)),
                                                  LGP.mean=round(mean(LGP,na.rm=T),1),
                                                  LGP.mode=getmode(LGP,na.rm=T),
                                                  LGP.median=round(median(LGP,na.rm=T),1),
                                                  LGP.min=suppressWarnings(min(LGP,na.rm=T)),
                                                  LGP.max=suppressWarnings(max(LGP,na.rm=T)),
                                                  LGP.sd=suppressWarnings(sd(LGP,na.rm=T)),
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
  
  
  LTAvg_SOS3<-Seasonal3[!is.na(Dekad.Season),list(Total.Seasons=.N,
                                                  SOS.mean=round(CircMean(m=SOS,interval=36,na.rm=T),1),
                                                  SOS.median=as.numeric(median(SOS,na.rm=T)),
                                                  SOS.mode=getmode(SOS,na.rm=T),
                                                  SOS.min=suppressWarnings(min(SOS,na.rm=T)),
                                                  SOS.max=suppressWarnings(max(SOS,na.rm=T)),
                                                  SOS.sd=suppressWarnings(sd(SOS,na.rm=T)),
                                                  EOS.mean=round(CircMean(m=EOS,interval=36,na.rm=T),1),
                                                  EOS.mode=getmode(EOS,na.rm=T),
                                                  EOS.median=round(median(EOS,na.rm=T),1),
                                                  EOS.min=suppressWarnings(min(EOS,na.rm=T)),
                                                  EOS.max=suppressWarnings(max(EOS,na.rm=T)),
                                                  EOS.sd=suppressWarnings(sd(EOS,na.rm=T)),
                                                  LGP.mean=round(mean(LGP,na.rm=T),1),
                                                  LGP.mode=getmode(LGP,na.rm=T),
                                                  LGP.median=round(median(LGP,na.rm=T),1),
                                                  LGP.min=suppressWarnings(min(LGP,na.rm=T)),
                                                  LGP.max=suppressWarnings(max(LGP,na.rm=T)),
                                                  LGP.sd=suppressWarnings(sd(LGP,na.rm=T)),
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
  
  # 8) Combine into a list and return data ####
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

SaveDir<-paste0(SOS_Dir,"/",FolderName)

if(!dir.exists(SaveDir)){
  dir.create(SaveDir,recursive=T)
}

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
