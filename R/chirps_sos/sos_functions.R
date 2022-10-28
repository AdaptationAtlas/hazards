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
#' `SOS_SeqMerge` takes a sequence created by `SOS_RSeason` where multiple non-NA values are split by sequences of NA values and in the corresponding aridity index (`AI`) data identifies the start and end point of non-NA values split by NA sequences of length `<=MaxGap`.
#'  Values from the start to end points are set to the first numeric value in the sequence. Values outwith the start and end points are set to NA.
#'  
#'  For example `NA NA NA NA NA  2  2  2  2  2 NA  3  3  3 NA NA  4 NA NA` with `MaxVal=1` becomes `NA NA NA NA NA  2  2  2  2  2  2  2  2  2 NA NA NA NA NA NA` and with `MaxVal=2` 
#'  becomes `NA NA NA NA NA  2  2  2  2  2  2  2  2  2  2  2  2 NA NA`.
#'
#' @param Seq  a vector of sequences split by NAs
#' @param AI a logical vector the same length as `Seq` indicating whether an aridity index threshold is `TRUE` or `FALSE`
#' @param MaxGap an integer value describing the maximum gap (number of NA values) allowed between non-NA values before the sequence breaks.
#' @param MinStartLen an integer value describing the minimum length of first sequence block, if the first sequence is too short and too separated from the next sequence it is removed. This is to remove false starts.
#' @param MaxStartSep an integer value describing the maximum separation of the first sequence block, if the first sequence is too short and too separated from the next sequence it is removed. This is to remove false starts.
#' @param ClipAI logical `T/F`, if `T` then `AI` values corresponding to the last non-NA value in `Seq` are all set to `F` and the sequence is halted by the last `F` value of AI.
#' @param S1.AI  a logical vector, when S1.AI==T then the corresponding AI value to the first value in a sequence of RAIN==T is set to `TRUE`. If `S1.AI` is set to `TRUE` in the functions that generate the `Seq` parameter it should also be set to `TRUE` here.
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
  conv=2*pi/interval
  x1 = Arg(mean(exp(conv*(m-1)*1i),na.rm=na.rm))
  x2 = x1/conv
  x3 = (x2 + interval) %% interval
  return(x3)
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
