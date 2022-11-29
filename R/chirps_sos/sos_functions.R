#' Find and Code Wettest Months
#'
#' `MaxMonth` takes 12 months of rainfall data and searches for two separate seasons with the highest rainfall.
#'
#' Rainfall data need to be a three month rolling average (n+2).
#' The `Pad` parameter sets the minimum separation distance between seasons in months (the maximum is 3).
#' The three months corresponding to the highest rainfall are coded as `1`, the next highest rainfall period
#' which does not overlap season 1 and is separated by `Pad` months is coded as `2`.
#' If two identical rainfall values exist in the 12 month series the earliest month will be chosen first.
#'
#' @param Rain a numeric vector (length = 12) of three month rolling average rainfall for a year (i.e 12 months).
#' @param Month  an integer vector (length = 12) of month numbers (1:12), if months are not in sequence the function should be able to handle this. If not supplied it is assumed months are in the order 1:12.
#' @param Pad an integer value (length = 1, values = 0:3) that sets the minimum separation distance between seasons in months (the maximum is 3).
#' @return `MaxMonth` returns a numeric vector (length = 12) with the highest rainfall season coded as `1` and the next highest encoded as `2`.
#' @export
SOS_MaxMonth<-function(Rain,Month,Pad){
  
  X<-rep(NA,12)
  
  N<-Month[Rain==max(Rain)]
  
  if(length(N)>1){N<-min(N)}
  
  N<-(N-Pad-2):(N+2+Pad)
  N[N>12]<-N[N>12]-12
  N[N<1]<-N[N<1]+12
  
  X[match(N,Month)]<-c(rep(99,Pad+2),rep(1,3),rep(99,Pad))
  
  Rain2<-Rain[is.na(X)]
  
  N2<-Month[is.na(X)][Rain2==max(Rain2)]
  if(length(N2)>1){N2<-min(N2)}
  N2<-N2:(N2+2)
  N2[N2>12]<-N2[N2>12]-12
  N2[N2<1]<-N2[N2<1]+12
  
  X[match(N2,Month)]<-rep(2,3)
  X[X==99]<-NA
  
  return(X)
  
}
#' Format date as yearly or monthly dekad (10-day period)
#'
#' `SOS_Dekad` transforms a class `Date` object to dekad of the year or month.
#'
#' This function is taken from Melanie Bacou's https://github.com/tidyverse/lubridate/issues/617#issuecomment-521393447
#'
#' @param x class `Date` vector to convert to dekad
#' @param type dekad of `month` (1:3) or `year` (1:36)
#' @inheritDotParams base::as.Date
#' @return integer dekad
#' @importFrom lubridate day
#' @examples
#' dekad(Sys.Date())
#' @export
SOS_Dekad <- function(Data,
                  type = c("month", "year"),
                  ...) {
  type <- match.arg(type)
  x <- as.Date(Data, ...)
  res <- ifelse(day(x) > 20,  3, ifelse(day(x) > 10, 2, 1))
  if(type == "year") res <- month(x)*3 + res - 3
  return(res)
}
#' Pad season sequences
#'
#' `SOS_SeasonPad` takes season sequences (i.e., values of 1 and 2 separated by NAs) and increases the sequence length in forward and backwards directions.
#'
#' @param Data an integer vector of seasons (sequences of `1` and `2` separated by NAs)
#' @param PadBack an integer value specifying number of places to increase sequence backwards. This should not be greater than the minimum separation of sequences.
#' @param PadForward an integer value specifying number of places to increase sequence forwards. This should not be greater than the minimum separation of sequences.
#' @return integer season
#' @export
SOS_SeasonPad<-function(Data,PadBack,PadForward){
  Data[is.na(Data)]<-99
  
  PadBack<-round(PadBack,0)
  PadForward<-round(PadForward,0)
  
  if(PadBack!=0){
    if(1 %in% Data){
      
    Seq1<-c(rep(99,PadBack),1)
    
    N <- which(Data == Seq1[1])
    N<-N[sapply(N, function(i) all(Data[i:(i+(length(Seq1)-1))] == Seq1))]
    N<-N[!is.na(N)]
    N<-c(N,rep(N,length(1:PadBack))+rep(1:PadBack,each=length(N)))
    N[N>length(Data)]<-NULL
    Data[N]<-1
    }
  
    if(2 %in% Data){
      Seq2<-c(rep(99,PadBack),2)
      
    N <- which(Data == Seq2[1])
    N<-N[sapply(N, function(i) all(Data[i:(i+(length(Seq2)-1))] == Seq2))]
    N<-N[!is.na(N)]
    N<-c(N,rep(N,length(1:PadBack))+rep(1:PadBack,each=length(N)))
    N[N>length(Data)]<-NULL
    Data[N]<-2
  }
  
  
  }
  
  if(PadForward!=0){
    if(1 %in% Data){
  Seq3<-c(1,rep(99,PadForward))
  
  N <- which(Data == Seq3[1])
  N<-N[sapply(N, function(i) all(Data[i:(i+(length(Seq3)-1))] == Seq3))]
  N<-N[!is.na(N)]
  N<-c(N,rep(N,length(1:PadForward))+rep(1:PadForward,each=length(N)))
  N[N>length(Data)]<-NULL
  Data[N]<-1
    }
  
  if(2 %in% Data){
    Seq4<-c(2,rep(99,PadForward))
    
    N <- which(Data == Seq4[1])
    N<-N[sapply(N, function(i) all(Data[i:(i+(length(Seq4)-1))] == Seq4))]
    N<-N[!is.na(N)]
    N<-c(N,rep(N,length(1:PadForward))+rep(1:PadForward,each=length(N)))
    N[N>length(Data)]<-NULL
    Data[N]<-2
  }
  }
  
  Data[Data==99]<-NA
  
  return(Data)
}
#' Give season sequences in a time-series unique values
#'
#' `SOS_UniqueSeq` takes season sequences (i.e., values of 1 and 2 separated by NAs) and gives each season a unique integer value.
#'
#' @param Data an integer vector of seasons (sequences of `1` and `2` separated by NAs)
#' @return integer unique sequence
#' @export
SOS_UniqueSeq<-function(Data){
  Y<-rle(Data)
  Seq<-rep(NA,length(Data))
  if(!all(is.na(Y$values))){
    Seq[!is.na(Data)]<-rep(1:sum(!is.na(Y$values)),Y$lengths[!is.na(Y$values)])
  }
  return(Seq)
}
#' Find rainy seasons in continuous dekadal sequences
#'
#' `SOS_RSeason` generates dekadal sequences where an initial rainfall condition (onset or SOS) has been met ending when an aridity condition is not met (end of season, EOS). It
#' takes two logical sequences of dekadal information `RAIN` and `AI` and determines rainy seasons from these.
#'
#' The `RAIN` parameter indicates where a starting condition of rainfall is `TRUE`. The `AI` parameter indicates an ending condition
#' where the aridity index is `FALSE`. The growing season occurs for sequences starting when `RAIN` is
#' `TRUE`. `AI` can be `T` or `F` for the starting dekad, but subsequent `AI` values must be `T` for the sequence to continue and the seasons ends when `AI` is `FALSE.`
#'
#' The typical WRSI onset date definition is a dekad with 25 mm followed by two dekads with a total of 20mm.
#' The typical WRSI end of season definition is when the long-term mean aridity index drops below 0.5 (i.e. mean potential evapotranspiration >= 2*mean rainfall).
#'
#' @param RAIN  a logical vector indicating where a starting condition of rainfall is `TRUE`
#' @param AI  a logical vector indicating whether an aridity index threshold is `TRUE` or `FALSE`
#' @param S1.AI  a logical vector, when S1.AI==T then the corresponding AI value to the first value in a sequence of RAIN==T is set to `TRUE`. If `S1.AI` is set to `TRUE` in the functions that generate the `Seq` parameter it should also be set to `TRUE` here.
#' @return integer sequence
#' @export
SOS_RSeason<-function(RAIN,AI,S1.AI){
  RAIN[is.na(RAIN)]<-F
  AI[is.na(AI)]<-F
  
  # Remove instances at start of AI==T sequence where RAIN==F
  SeqFun2<-function(RAIN,SEQ){
    N<-which(RAIN==T)[1] # Find first instance where both are true, if index is >1 then this means that there is an instance of F/T
    SEQ<-rep(SEQ,length(RAIN))
    
    # If all instance of AI & RAIN are F then the rain threshold was not met and there is no sequence.
    if(is.na(N)){
      SEQ<-NA
    }else{
      if(N>1){
        SEQ[1:(N-1)]<-NA  # Set T/F values to NA
      }
    }
    return(SEQ)
  }
  
  # If S1.AI==T then the corresponding AI value to the first value in a sequence of RAIN==T is also set to TRUE
  if(S1.AI){
    RainSeq<-rle(RAIN)
    CSrain<-cumsum(RainSeq$length)
    CSrain[1]<-1
    if(length(CSrain)>1){
    CSrain[2:length(CSrain)]<-CSrain[1:(length(CSrain)-1)]+1
    }
    
    N<-which(RainSeq$values==T)
    
    AI[CSrain[N]]<-T
  }
  
  
  # Find T sequences
  SeqLen<-rle(AI)
  Lengths<-SeqLen$lengths
  Values<-SeqLen$values
  
  # If all values = F then there was no sequence
  if(all(Values==F)){
    Seq<-as.integer(rep(NA,length(RAIN)))
  }else{
    Y<-rep(NA,length(AI))
    Y[AI==T & !is.na(AI)]<-rep(1:sum(Values,na.rm = T),Lengths[Values==T & !is.na(Values)])
    Y[1]<-T
    Seq<-data.table(AI=AI,RAIN=RAIN,SEQ=Y)
    
    Seq<-Seq[!is.na(SEQ),SEQ:=SeqFun2(RAIN,SEQ),by=SEQ][,SEQ]
  }
  
  return(Seq)
}
#' Merge sequences for a continuous series of dekadal data
#'
#' `SOS_SeqMerge` takes a sequence created by `SOS_RSeason` where multiple non-NA values are split by sequences of NA values and 
#' in the corresponding aridity index (`AI`) data identifies the start and end point of non-NA values split by NA sequences of length `<=MaxGap`.
#'  Values from the start to end points are set to the first numeric value in the sequence. Values outwith the start and end points are set to NA.
#'  
#'  For example `NA NA NA NA NA  2  2  2  2  2 NA  3  3  3 NA NA  4 NA NA` with `MaxVal=1` becomes 
#'  `NA NA NA NA NA  2  2  2  2  2  2  2  2  2 NA NA NA NA NA NA` and with `MaxVal=2` 
#'  becomes `NA NA NA NA NA  2  2  2  2  2  2  2  2  2  2  2  2 NA NA`.
#'
#' @param Seq  a vector of sequences split by NAs
#' @param AI a logical vector the same length as `Seq` indicating whether an aridity index threshold is `TRUE` or `FALSE`
#' @param MaxGap an integer value describing the maximum gap (number of NA values) allowed between non-NA values before the sequence breaks.
#' @param MinStartLen an integer value describing the minimum length of first sequence block, if the first sequence is too short and too 
#' separated from the next sequence it is removed. This is to remove false starts.
#' @param MaxStartSep an integer value describing the maximum separation of the first sequence block, if the first sequence is too short and too
#'  separated from the next sequence it is removed. This is to remove false starts.
#' @param ClipAI logical `T/F`, if `T` then `AI` values corresponding to the last non-NA value in `Seq` are all set to `F` and the sequence is 
#' halted by the last `F` value of AI.
#' @param S1.AI  a logical vector, when S1.AI==T then the corresponding AI value to the first value in a sequence of RAIN==T is set to `TRUE`. 
#' If `S1.AI` is set to `TRUE` in the functions that generate the `Seq` parameter it should also be set to `TRUE` here.
#' @return a vector where sequences are merged across NAs and set to the first value of the series. Leading and trailing NAs remain.
#' @export
SOS_SeqMerge<-function(Seq,AI,MaxGap,MinStartLen,MaxStartSep,ClipAI,S1.AI){
  
  # If only one sequence is present then nothing needs to happen
  if(!(all(is.na(Seq))|length(unique(Seq))==1)){
    
    # If S1.AI==T then set the corresponding AI value to the first value in a sequence of RAIN==T to TRUE
    if(S1.AI){
      Seq.Seq<-rle(Seq)
      CSSeq<-cumsum(Seq.Seq$length)
      CSSeq[1]<-1
      CSSeq[2:length(CSSeq)]<-CSSeq[1:(length(CSSeq)-1)]+1
      
      N<-which(!is.na(Seq.Seq$values))
      
      AI[CSSeq[N]]<-T
    }
    
    # If starting separation cannot be less than the max allowable gap
    if(MaxGap>MaxStartSep){
      MaxStartSep<-MaxGap
    }
    
    # Set initial AI values where season start is false to false
    N2<-which(!is.na(Seq))
    if((N2[1]-1)>0){
      AI[1:(N2[1]-1)]<-F
    }
    
    NSeq<-length(unique(Seq[!is.na(Seq)]))
    if(!ClipAI & NSeq==1){
      N1<-which(!is.na(Seq))[1]
      N1<-c(rep(F,N1-1),rep(T,length(Seq)-N1+1))
      Seq[AI==T & N1==T]<-Seq[!is.na(Seq)][1]
    }
    
    Seq2<-Seq
    Seq2[is.na(Seq2)]<-99999
    X<-rle(Seq2)
    SVal<-X$values
    SLen<-X$lengths
    SN<-which(SVal!=99999)
    
    if(length(SN)>1){
      
      if(N2[length(N2)]<length(AI) & ClipAI){
        AI[(N2[length(N2)]+1):length(AI)]<-F
      }
      
      X<-rle(AI)
      Val<-X$values
      Length<-X$lengths
      Z<-rep(NA,length(AI))
      N<-which(Val==T)
      
      A<-N[1]:N[length(N)]
      # Find NA sequences in between non-NA seqs
      B<-A[!A %in% N]
      
      # Remove false starts
      if(Length[N[1]]<=MinStartLen & Length[B[1]]>MaxStartSep & NSeq>1){
        i<-which(AI==T)[1]
        
        # When is next rainfall event?
        j<-cumsum(SLen)[SN[2]-1]
        AI[i:j]<-F
        
        X<-rle(AI)
        Val<-X$values
        Length<-X$lengths
        Z<-rep(NA,length(AI))
        N<-which(Val==T)
        A<-N[1]:N[length(N)]
        # Find NA sequences in between non-NA seqs
        B<-A[!A %in% N]
      }
      
      if(length(N)>1){
        # How many consecutive B sequences <= MaxGap (starting from the first within-seq NA break)
        C<-Length[B]<=MaxGap
        if(C[1]==T){
          i<-rle(C)$lengths[1]+1
        }else{
          i<-1
        }
        
        j<-if(N[1]!=1){Length[1]+1}else{1} # Start of non-NA sequence
        k<-sum(Length[1:N[i]]) # End of non-NA sequence
        
        Z[j:k]<-Seq[!is.na(Seq)][1]
        
      }else{
        Z[AI==T]<-Seq[!is.na(Seq)][1]
      }
      
      return(Z)
      
    }else{
      
      return(Seq) 
    }
  }else{
    return(Seq)
  }
}
#' Merge sequences for a long term average sequence of dekadal data
#' @export
LTSeqMerge<-function(Seq){
  if(length(unique(Seq[!is.na(Seq)]))>1){
    Seq[is.na(Seq)]<-"X"
    
    XX<-Seq[1]=="X" & tail(Seq,1)=="X"
    if(XX){
      Seq2<-rle(Seq)
      Adj<-1:36-Seq2$lengths[1]
      Adj[Adj<1]<-Adj[Adj<1]+36
      
      AdjBack<-1:36+Seq2$lengths[1]
      AdjBack[AdjBack>36]<-AdjBack[AdjBack>36]-36
      
      SeqAdj<-Seq[Adj]
      Seq<-rle(SeqAdj)
    }else{
      Seq<-rle(Seq)
    }
    X<-which(Seq$values=="X")
    N1<-X[which(Seq$lengths[X]==1)]
    if(length(N1)>0){
      N2<-N1-1
      N2[N2<1]<-length(Seq$values)
      Seq$values[N1]<-Seq$values[N2]
      Seq<-rep(Seq$values,Seq$lengths)
      
      if(XX){
        Seq<-Seq[AdjBack]
      }
      
      Seq[Seq!="X"]<-1
      Seq<-rle(Seq)
      Seq$values[Seq$values!="X"]<-1:length(Seq$values[Seq$values!="X"])
      
      Seq<-suppressWarnings(as.numeric(rep(Seq$values,Seq$lengths)))
      
    }else{
      Seq<-suppressWarnings(as.numeric(rep(Seq$values,Seq$lengths)))
      
      if(XX){
        Seq<-Seq[AdjBack]
      }
    }
    
    Seq[Seq=="X"]<-NA
    
    if(!is.na(Seq[1]) & !is.na(Seq[length(Seq)])){
      Seq<-rle(Seq)
      Seq$values[length(Seq$values)]<-Seq$values[1]
      Seq<-rep(Seq$values,Seq$lengths)
    }
  }
  return(Seq)
}
# Find rainy seasons in long-term average dekadal sequence
#' @export
SOS_RSeasonLT<-function(RAIN,AI,S1.AI){
  
  # Remove instances at start of AI==T sequence where RAIN==F
  SeqFun2<-function(RAIN,SEQ){
    RAIN[is.na(RAIN)]<-F
    N<-which(RAIN==T)[1] # Find first instance where both are true, if index is >1 then this means that there is an instance of F/T
    SEQ<-rep(SEQ,length(RAIN))
    
    # If all instance of AI & RAIN are F then the rain threshold was not met and there is no sequence.
    if(is.na(N)){
      SEQ<-NA
    }else{
      if(N>1){
        SEQ[1:(N-1)]<-NA  # Set T/F values to NA
      }
    }
    return(SEQ)
  }
  
  # If start and end of sequence are both FALSE 
  XX<-RAIN[1]==F & tail(RAIN,1)==F
  
  if(XX){
    Seq2<-rle(RAIN)
    Adj<-1:36-Seq2$lengths[1]
    Adj[Adj<1]<-Adj[Adj<1]+36
    
    AdjBack<-1:36+Seq2$lengths[1]
    AdjBack[AdjBack>36]<-AdjBack[AdjBack>36]-36
    
    RAIN2<-RAIN[Adj]
  }else{
    RAIN2<-RAIN
  }
  
  if(S1.AI){
  RainSeq<-rle(RAIN2)
  CSrain<-cumsum(RainSeq$length)
  CSrain[1]<-1
  CSrain[2:length(CSrain)]<-CSrain[1:(length(CSrain)-1)]+1
  
  N<-which(RainSeq$values==T)

  N3<-CSrain[N]
  
  if(XX){
    N3<-Adj[N3]
  }
  
  AI[N3]<-T
  }
  
  # Find T sequences
  SeqLen<-rle(AI)
  Lengths<-SeqLen$lengths
  Values<-SeqLen$values
  
  # If all values = F then there was no sequence
  if(all(Values==F)){
    Seq<-as.integer(rep(NA,length(RAIN)))
  }else{
    Y<-rep(NA,length(AI))
    Y[AI==T & !is.na(AI)]<-rep(1:sum(Values,na.rm = T),Lengths[Values==T & !is.na(Values)])
    Y[1]<-T
    Seq<-data.table(AI=AI,RAIN=RAIN,SEQ=Y)
    
    Seq<-Seq[!is.na(SEQ),SEQ:=SeqFun2(RAIN,SEQ),by=SEQ][,SEQ]
  }
  
  if(!is.na(Seq[1]) & !is.na(Seq[length(Seq)])){
    Seq<-rle(Seq)
    Seq$values[length(Seq$values)]<-Seq$values[1]
    Seq<-rep(Seq$values,Seq$lengths)
  }
  
  return(Seq)
}
# Order a dekadal sequence that crosses the year end boundary
#' @export
OrderDekadSeq<-function(x){
  x<-sort(x)
  xo<-x[1]:(x[1]+length(x)-1)
  if(any(xo!=x)){
    x<-c(x[xo!=x],x[xo==x])
  }
  return(x)
}
# Calculate circular mean of times
#' @export
CircMean <- function(m,interval,na.rm=T){
  conv <- 2*pi/36
  m <- mean(exp(conv*(m[!is.na(m)])*1i))
  m <- interval+Arg(m)/conv%%interval  ## 'direction', i.e. average month
  m <- (m + interval) %% interval
  return(m)
}
# Calculate mean deviance from average value of circular vector of values
#' @export
mean_deviance<-function(data,interval=36,na.rm=T,fun="mean",stat="mean"){
  data<-data[!is.na(data)]
  if(length(data)>0){
    a<-if(fun=="mean"){
      CircMean(data,interval,na.rm)
    }else{
      getmode(data,na.rm)
    }
    b<-sapply(data,FUN=function(x){CicularDist(Val1=x,Val2=a,interval=interval)})
    x<-if(stat=="mean"){
      mean(b)
    }else{
      if(stat=="sd"){
        sd(b)
      }}
    return(x)
  }else{
    return(NA)
  }
  
}

# Applies a function to subsets of a given dataset. https://gist.github.com/stillmatic/fadfd3269b900e1fd7ee
#' @export
slide_apply <- function (data, window, step = 1, fun) 
{
  fun <- match.fun(fun)
  total <- length(data)
  window <- abs(window)
  spots <- seq(from = 1, to = (total - window + 1), by = abs(step))
  result <- rep(NA, length(spots))
  for (i in 1:length(spots)) {
    result[window + i - 1] <- fun(data[spots[i]:(spots[i] + 
                                                   window - 1)])
  }
  
  result<-result[c(which(!is.na(result)),which(is.na(result)))]-data
  return(result)
}
# Applies a function to subsets of a given dataset. https://gist.github.com/stillmatic/fadfd3269b900e1fd7ee
#' @export
slide_apply2 <- function (data, window, step = 1, fun) 
{
  fun <- match.fun(fun)
  total <- length(data)
  window <- abs(window)
  spots <- seq(from = 1, to = (total - window + 1), by = abs(step))
  result <- rep(NA, length(spots))
  for (i in 1:length(spots)) {
    result[window + i - 1] <- fun(data[spots[i]:(spots[i] + 
                                                   window - 1)])
  }
  
  result<-result[c(which(!is.na(result)),which(is.na(result)))]
  return(result)
}
# This function takes two angles and calculates the difference between then as such we need to convert dekads 1:36 to degrees 1:360 (and back again)
#' @export
CicularDist<-function(Val1,Val2,interval){
  Max<-max(Val1,Val2)
  Min<-min(Val1,Val2)
  
  m<-min(Max-Min,Min-Max+interval)
  
  return(m)
}
# Calculate mode
#' @export
getmode <- function(v,na.rm) {
  if(na.rm){
    v<-v[!is.na(v)]
  }
uniqv <- unique(v)
uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Dekad to month
#' @export
dkd2mth<-function(x){
  dkd<-1:36
  mth<-rep(1:12,each=3)
  return(mth[match(x,dkd)])
}
# Read raster into memory
#' https://rdrr.io/github/rspatial/terra/src/R/read.R
#' @export

readAll <- function(x) {
  ok <- x@ptr$readAll()
  x <- messages(x, "readAll")
  invisible(ok)
}

setMethod("readStart", signature(x="SpatRaster"),
          function(x) {
            success <- x@ptr$readStart()
            messages(x, "readStart")
            if (!success) error("readStart,SpatRaster", "cannot open file for reading")
            invisible(success)
          }
)

setMethod("readStart", signature(x="SpatRasterDataset"),
          function(x) {
            success <- x@ptr$readStart()
            messages(x, "readStart")
            if (!success) error("readStart,SpatRasterDataset", "cannot open file for reading")
            invisible(success)
          }
)

setMethod("readStop", signature(x="SpatRaster"),
          function(x) {
            success <- x@ptr$readStop()
            messages(x, "readStop")
            invisible(success)
          }
)

setMethod("readStop", signature(x="SpatRasterDataset"),
          function(x) {
            success <- x@ptr$readStop()
            messages(x, "readStop")
            invisible(success)
          }
)
# Terra error messages - error
#' https://rdrr.io/github/rspatial/terra/src/R/messages.R
#' @export
error <- function(f, emsg="", ...) {
  stop("[", f, "] ", emsg, ..., call.=FALSE)
}
# Terra error messages - warn
#' https://rdrr.io/github/rspatial/terra/src/R/messages.R
#' @export
warn <- function(f, wmsg="", ...) {
  warning("[", f, "] ", wmsg, ..., call.=FALSE)
}
# Terra error messages - messages
#' https://rdrr.io/github/rspatial/terra/src/R/messages.R
#' @export
messages <- function(x, f="") {
  #g <- gc(verbose=FALSE)
  if (methods::.hasSlot(x, "ptr")) {
    if (x@ptr$has_warning()) {
      warn(f, paste(unique(x@ptr$getWarnings()), collapse="\n"))
    }
    if (x@ptr$has_error()) {
      error(f, x@ptr$getError())
    }
  } else {
    if (x$has_warning()) {
      warn(f, paste(unique(x$getWarnings()), collapse="\n"))
    }
    if (x$has_error()) {
      error(f, x$getError())
    }
  }
  x
}
# SOS analysis function
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
#' @export
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
    DATA<-DATA[,list(Rain.Dekad=sum(Rain),ETo=sum(ETo),AI=mean(AI)),by=list(Index,Year,Dekad,Rain.Season) # Sum rainfall and take mean aridity, Index by dekad  (within year and rain season)
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
  ]
  
  return(DATA)
  
}
# Function to add a similarity field (proportion of SOS dekads which are the same as the most frequent SOS dekads). Used in SOS_Wrap function.
#' @export
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
# Function to calculate season separation. Used in SOS_Wrap function.
#' @export
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
# Create a wrapper for data.table operations
#' @param Season2.Prop
#' @param MinLength
#' @param RollBack
#' @export
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
  
  CLIM.Dekad[!(is.na(AI.Seq)|is.na(Dekad.Seq)),Start.Year:=Year[1],by=list(Index,Dekad.Seq) # Add starting year for seasons
             ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq)  # Add total rainfall for season # Add total rainfall for season
               ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.ETo:=sum(ETo),by=list(Index,Dekad.Seq,AI.Seq)] # Add total ETo for season # Add total rainfall for season
  
  # 2) Calculate Seasonal Values ####
  Len<-CLIM.Dekad[,length(unique(Year))]
  Seasonal<-unique(CLIM.Dekad[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,SOS,EOS,LGP,Dekad.Season,Tot.Rain)])
  Seasonal[!is.na(Dekad.Season),Seasons.Count:=.N,by=list(Index,Dekad.Season)
  ][,Season2Prop:=Seasons.Count/Len,by=Index]
  
  Seasonal[,Seasons:=length(unique(Dekad.Season)),by=Index]
  
  # 3) Roll back SOS where SOS is fixed ####
  if(RollBack==T){

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
      
      CLIM.Dekad2[!(is.na(AI.Seq)|is.na(Dekad.Seq)),Start.Year:=Year[1],by=list(Index,Dekad.Seq) # Add starting year for seasons
      ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq)  # Add total rainfall for season # Add total rainfall for season
      ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.ETo:=sum(ETo),by=list(Index,Dekad.Seq,AI.Seq)] # Add total ETo for season # Add total rainfall for season
      
      
      CLIM.Dekad<-rbind(CLIM.Dekad[!Index %in% Sites],CLIM.Dekad2)
    }
  }
  # 3.2) Scenario 2: SOS fixed and two seasons present #####
  if(RollBack==T){
    # Subset data
    Data<-CLIM.Dekad[Index %in% X[Seasons==2,Index]]
    
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
      
      CLIM.Dekad1[!(is.na(AI.Seq)|is.na(Dekad.Seq)),Start.Year:=Year[1],by=list(Index,Dekad.Seq) # Add starting year for seasons
      ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq)  # Add total rainfall for season # Add total rainfall for season
      ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.ETo:=sum(ETo),by=list(Index,Dekad.Seq,AI.Seq)] # Add total ETo for season # Add total rainfall for season
      
      
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
      CLIM.Dekad1<-DATA[Index %in% Sites,list(Rain.Dekad=sum(Rain),ETo=sum(ETo),AI=mean(AI)),by=list(Index,Year,Dekad,Rain.Season)]
      
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
      
      CLIM.Dekad1[!(is.na(AI.Seq)|is.na(Dekad.Seq)),Start.Year:=Year[1],by=list(Index,Dekad.Seq) # Add starting year for seasons
      ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq)  # Add total rainfall for season # Add total rainfall for season
      ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.ETo:=sum(ETo),by=list(Index,Dekad.Seq,AI.Seq)] # Add total ETo for season # Add total rainfall for season
      
      
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
        ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq) # Add total rainfall for season
        ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.ETo:=sum(ETo),by=list(Index,Dekad.Seq,AI.Seq)] # Add ETo rainfall for season
      }
      
      CLIM.Dekad<-rbind(CLIM.Dekad[!Index %in% Sites],CLIM.Dekad1)
      Clim.Dekad1<-NULL
    }
  }
  # 4.3) Calculate seasonal values #####
  
  Seasonal2<-unique(CLIM.Dekad[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,SOS,EOS,LGP,Dekad.Season,Tot.Rain,Tot.ETo)])
  # Remove second seasons that are too short
  if(!is.na(MinLength)){
    Seasonal2<-Seasonal2[!(Dekad.Season==2 & LGP<MinLength)]
  }
  Seasonal2<-Seasonal2[!is.na(Dekad.Season),Seasons.Count:=.N,by=list(Index,Dekad.Season)][,Season2Prop:=Seasons.Count/Len]
  
  # Remove second seasons that are present for less than 1/3 the time of first seasons
  if(!is.na(Season2.Prop)){
    Seasonal2<-Seasonal2[Dekad.Season==2 & Season2Prop>Season2.Prop]
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
    ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.Rain:=sum(Rain.Dekad),by=list(Index,Dekad.Seq,AI.Seq) # Add total rainfall for season
    ][!(is.na(AI.Seq)|is.na(Dekad.Seq)),Tot.ETo:=sum(ETo),by=list(Index,Dekad.Seq,AI.Seq)] # Add ETo rainfall for season
    
    Seasonal3<-unique(CLIM.Dekad1[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,SOS,EOS,LGP,Dekad.Season,Tot.Rain,Tot.ETo)])
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
      Seasonal3<-unique(CLIM.Dekad.3[!(is.na(Dekad.Season)|is.na(Start.Year)),list(Index,Start.Year,Dekad.Season,SOS,EOS,LGP,Tot.Rain,Tot.ETo)])
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
                                                  SOS.mean.dev=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mean",stat="mean"),2),
                                                  SOS.mean.dev.sd=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mean",stat="sd"),2),
                                                  SOS.mode=getmode(SOS,na.rm=T),
                                                  SOS.mode.dev=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mode",stat="mean"),2),
                                                  SOS.mode.dev.sd=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mode",stat="sd"),2),
                                                  SOS.min=suppressWarnings(min(SOS,na.rm=T)),
                                                  SOS.max=suppressWarnings(max(SOS,na.rm=T)),
                                                  EOS.mean=round(CircMean(m=EOS,interval=36,na.rm=T),1),
                                                  EOS.mean.dev=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mean",stat="mean"),2),
                                                  EOS.mean.dev.sd=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mean",stat="sd"),2),
                                                  EOS.mode=getmode(EOS,na.rm=T),
                                                  EOS.mode.dev=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mode",stat="mean"),2),
                                                  EOS.mode.dev.sd=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mode",stat="sd"),2),
                                                  EOS.min=suppressWarnings(min(as.integer(EOS),na.rm=T)),
                                                  EOS.max=suppressWarnings(max(as.integer(EOS),na.rm=T)),
                                                  LGP.mean=round(mean(LGP,na.rm=T),1),
                                                  LGP.mode=getmode(LGP,na.rm=T),
                                                  LGP.min=suppressWarnings(min(LGP,na.rm=T)),
                                                  LGP.max=suppressWarnings(max(LGP,na.rm=T)),
                                                  LGP.sd=round(suppressWarnings(sd(LGP,na.rm=T)),2),
                                                  Tot.Rain.mean=round(mean(Tot.Rain,na.rm=T),1),
                                                  Tot.Rain.sd=round(sd(Tot.Rain,na.rm=T),2),
                                                  Tot.ETo.mean=round(mean(Tot.ETo,na.rm=T),1),
                                                  Tot.ETo.sd=round(sd(Tot.ETo,na.rm=T),2),
                                                  Balance.mean=round(mean(Tot.Rain-Tot.ETo,na.rm=T),1),
                                                  Balance.sd=round(sd(Tot.Rain-Tot.ETo,na.rm=T),2),
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
                                                  SOS.mean.dev=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mean",stat="mean"),2),
                                                  SOS.mean.dev.sd=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mean",stat="sd"),2),
                                                  SOS.mode=getmode(SOS,na.rm=T),
                                                  SOS.mode.dev=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mode",stat="mean"),2),
                                                  SOS.mode.dev.sd=round(mean_deviance(data=SOS,interval=36,na.rm=T,fun="mode",stat="sd"),2),
                                                  SOS.min=suppressWarnings(min(SOS,na.rm=T)),
                                                  SOS.max=suppressWarnings(max(SOS,na.rm=T)),
                                                  EOS.mean=round(CircMean(m=EOS,interval=36,na.rm=T),1),
                                                  EOS.mean.dev=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mean",stat="mean"),2),
                                                  EOS.mean.dev.sd=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mean",stat="sd"),2),
                                                  EOS.mode=getmode(EOS,na.rm=T),
                                                  EOS.mode.dev=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mode",stat="mean"),2),
                                                  EOS.mode.dev.sd=round(mean_deviance(data=EOS,interval=36,na.rm=T,fun="mode",stat="sd"),2),
                                                  EOS.min=suppressWarnings(min(EOS,na.rm=T)),
                                                  EOS.max=suppressWarnings(max(EOS,na.rm=T)),
                                                  LGP.mean=round(mean(LGP,na.rm=T),1),
                                                  LGP.mode=getmode(LGP,na.rm=T),
                                                  LGP.min=suppressWarnings(min(LGP,na.rm=T)),
                                                  LGP.max=suppressWarnings(max(LGP,na.rm=T)),
                                                  LGP.sd=round(suppressWarnings(sd(LGP,na.rm=T)),2),
                                                  Tot.Rain.mean=round(mean(Tot.Rain,na.rm=T),1),
                                                  Tot.Rain.sd=round(sd(Tot.Rain,na.rm=T),1),
                                                  Tot.ETo.mean=round(mean(Tot.ETo,na.rm=T),1),
                                                  Tot.ETo.sd=round(sd(Tot.ETo,na.rm=T),2),
                                                  Balance.mean=round(mean(Tot.Rain-Tot.ETo,na.rm=T),1),
                                                  Balance.sd=round(sd(Tot.Rain-Tot.ETo,na.rm=T),2),
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

