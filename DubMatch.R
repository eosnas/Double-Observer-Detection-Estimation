#' A function to determine matches from a double observer trail
#'
#' @param front.data A data frame or character path to a readable CSV file for front observer.
#' @param rear.data As above but for rear observer. 
#' @param combined.data A data frame or path to both front and rear observers. 
#' @param front.obs A character string that identifies the front observer. 
#' @param rear.obs A character string that identifies the rear observer. 
#' @param tdist Numeric scalar that defines the upper time distance between matched observations. 
#' @param sdist Numeric scalar that defines the upper spatial distance between matched observations. 
#' @param lag Time lag to apply to all observations of the rear observer, see details.  
#' @param open The threshold in group size below which exact matches are required. 
#' 
#' @return
#' @export
#'
#' @examples
DubMatch = function(front.data=NA, rear.data=NA, combined.data=NA, front.obs, rear.obs, 
                    tdist=5, sdist=250, lag=0, open=5){
  
  #check if combined data was given, if not, read in each file
  
  if(class(front.data)=="data.frame"){
    
    f=front.data[front.data$obs==front.obs, ]
    
  }
  
  if(class(rear.data)=="data.frame"){
    
    r=rear.data[rear.data$obs==rear.obs, ]
    
    
  }
  
  if(class(front.data)=="character"){
    
    data=read.csv(front.data, header=TRUE, stringsAsFactors = FALSE)
    f=data[data$obs==front.obs, ]
    
    
  }
  
  if(class(rear.data)=="character"){
    
    data=read.csv(rear.data, header=TRUE, stringsAsFactors = FALSE)
    r=data[data$obs==rear.obs, ]
    
    
  }
  
  if(class(combined.data)=="character"){
    
    data=read.csv(combined.data, header=TRUE, stringsAsFactors = FALSE)
    f=data[data$obs==front.obs, ]
    r=data[data$obs==rear.obs, ]
    
  }
  
  if(class(combined.data)=="data.frame"){
    
    f=combined.data[combined.data$obs==front.obs, ]
    r=combined.data[combined.data$obs==rear.obs, ]
    
    
  }
  
  f$matched=rep(0, dim(f)[1])
  r$matched=rep(0, dim(r)[1])
  
  
  f.tran=unique(f$tran)
  r.tran=unique(r$tran)
  
  common=as.character(f.tran[!is.na(match(f.tran, r.tran))])
  
  f=f[f$tran %in% common,]
  r=r[r$tran %in% common,]
  
  #empty data frame to populate with matches
  
  matches=data.frame(yr=integer(),
                     tran=character(),
                     ch=character(),
                     sppn=character(),
                     grp=integer(),
                     unit=character(),
                     front=character(),
                     rear=character(),
                     crew=character(), stringsAsFactors = FALSE
  )
  
  
  for (i in 1:length(f$yr)){
    
    s <- findNear(st=as.numeric(f[i,c("long", "lat", "ctime", "da", "mo", "yr")]), 
                  y=r[,c("long", "lat", "ctime", "da", "mo", "yr")], 
                  sdist=sdist, tdist=tdist, tlag=lag)
    pick <- s & r$matched!= 1 & r$tran == f$tran[i] & r$sppn == f$sppn[i] #& r$unit==f$unit[i]
    m <- r[pick, ]
  
    if ( length(m$yr) == 0 ){
      
      newline=data.frame(yr=f$yr[i],
                         tran=f$tran[i],
                         ch="10",
                         sppn=f$sppn[i],
                         grp=f$grp[i],
                         unit=f$unit[i],
                         front=f$obs[i],
                         rear=r$obs[1],
                         crew=paste(f$obs[i], r$obs[1], sep=""), stringsAsFactors = FALSE
      )
      
      matches=rbind(matches, newline)
      f$matched[i] = 1
    }
    
    if ( length(m$yr) == 1 ){
      if( m$unit==f$unit[i] ){
        newline=data.frame(yr=f$yr[i],
                           tran=f$tran[i],
                           ch="11",
                           sppn=f$sppn[i],
                           grp= (f$grp[i] + m$grp[1])/2,
                           unit=f$unit[i],
                           front=f$obs[i],
                           rear=m$obs[1],
                           crew=paste(f$obs[i], m$obs[1], sep=""), stringsAsFactors = FALSE
        )
        
        matches=rbind(matches, newline)
        f$matched[i] = 1
        r$matched[which(pick == TRUE)] = 1
      } else{
        if( m$unit=="open" & f$unit=="single" & m$grp==f$grp ){
          newline=data.frame(yr=f$yr[i],
                             tran=f$tran[i],
                             ch="11",
                             sppn=f$sppn[i],
                             grp= (f$grp[i] + m$grp[1])/2,
                             unit=f$unit[i],
                             front=f$obs[i],
                             rear=m$obs[1],
                             crew=paste(f$obs[i], m$obs[1], sep=""), stringsAsFactors = FALSE
          )
          
          matches=rbind(matches, newline)
          f$matched[i] = 1
          r$matched[which(pick == TRUE)] = 1
        }
        if( m$unit=="single" & f$unit=="open" & m$grp==f$grp ){
          newline=data.frame(yr=f$yr[i],
                             tran=f$tran[i],
                             ch="11",
                             sppn=f$sppn[i],
                             grp= (f$grp[i] + m$grp[1])/2,
                             unit=f$unit[i],
                             front=f$obs[i],
                             rear=m$obs[1],
                             crew=paste(f$obs[i], m$obs[1], sep=""), stringsAsFactors = FALSE
          )
          
          matches=rbind(matches, newline)
          f$matched[i] = 1
          r$matched[which(pick == TRUE)] = 1
        }
      }
    }
    #need to add open/single comparison here for length(m) > 1
    if ( length(m$yr) > 1 ){
      if( !is.na(f$long[1]) ){
        sdist2 = sapply(1:(dim(m)[1]), function(x){
          hav.dist(long1=f[i,"long"], lat1=f[i,"lat"], long2=m[x,"long"], lat2=m[x,"lat"])}) #in meters
        tdist2 = sapply(1:(dim(m)[1]), function(x){f[i,"ctime"] - m[x, "ctime"]}) #in seconds
        wdist = sqrt( (sdist2/45)^2 + tdist2^2 ) #assume moving at 45 meters per second
      } else{
        wdist = sapply(1:(dim(m)[1]), function(x){f[i,"ctime"] - m[x, "ctime"]}) #in seconds
      }
      mpick = which(wdist==min(wdist))[1]
      newline=data.frame(yr=f$yr[i],
                         tran=f$tran[i],
                         ch="11",
                         sppn=f$sppn[i],
                         grp= (f$grp[i] + m$grp[mpick])/2,
                         unit=f$unit[i],
                         front=f$obs[i],
                         rear=m$obs[mpick],
                         crew=paste(f$obs[i], m$obs[mpick], sep=""), stringsAsFactors = FALSE
      )
      matches=rbind(matches, newline)
      f$matched[i] = 1
      r$matched[which(pick == TRUE)][mpick] = 1
    }
    
  }
  

  newline=data.frame(yr=r$yr[r$matched!=1],
                     tran=r$tran[r$matched!=1],
                     ch=rep("01", sum(r$matched!=1)),
                     sppn=r$sppn[r$matched!=1],
                     grp=r$grp[r$matched!=1],
                     unit=r$unit[r$matched!=1],
                     front=rep(f$obs[1], sum(r$matched!=1)),
                     rear=r$obs[r$matched!=1],
                     crew=paste(rep(f$obs[1], sum(r$matched!=1)), r$obs[r$matched!=1], sep=""), stringsAsFactors = FALSE
  )
  
  matches=rbind(matches, newline)
  
  return(matches)
  
}