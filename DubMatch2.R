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
#' 
#' @return
#' @export
#'
#' @examples
DubMatch = function(front.data=NA, rear.data=NA, combined.data=NA, front.obs, rear.obs, 
                    tdist=5, sdist=250, lag=0){
  #Notes:
  #expand lines with >1 single to multiple lines?
  #add match of 2 single to 2 flkdrake
  #need to add matching of flkdrakes with opens? 

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
  
  
  # f.tran=unique(f$tran)
  # r.tran=unique(r$tran)
  # 
  # common=as.character(f.tran[!is.na(match(f.tran, r.tran))])
  # 
  # f=f[f$tran %in% common,]
  # r=r[r$tran %in% common,]
  
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
  
  
  for (i in 1:length(f$yr) ){
    
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
    
    if ( length(m$yr) > 0 ){
      #pick based on unit
      unitPick <- m$unit == f$unit[i]
      if( sum(unitPick) > 0 ){
        #pick based on grp size
        grpPick <- m$grp[unitPick] == f$grp[i]
        if( sum(grpPick) > 0 ){
          if( !is.na(f$long[1]) ){
            sdist2 = sapply(1:(dim(m[unitPick,][grpPick,])[1]), function(x){
              hav.dist(long1=f[i,"long"], lat1=f[i,"lat"], long2=m[unitPick,][grpPick,][x,"long"], 
                       lat2=m[unitPick,][grpPick,][x,"lat"])}) #in meters
            tdist2 = sapply(1:(dim(m[unitPick,][grpPick,])[1]), function(x){f[i,"ctime"] - m[unitPick,][grpPick,][x, "ctime"]}) #in seconds
            wdist = sqrt( (sdist2/45)^2 + tdist2^2 ) #assume moving at 45 meters per second
          } else{
            wdist = sapply(1:(dim(m[unitPick,][grpPick,])[1]), function(x){f[i,"ctime"] - m[unitPick,][grpPick,][x, "ctime"]}) #in seconds
          }
          mpick = which(wdist==min(wdist))[1]
          newline=data.frame(yr=f$yr[i],
                             tran=f$tran[i],
                             ch="11",
                             sppn=f$sppn[i],
                             grp= f$grp[i], #grp size is equal
                             unit=f$unit[i],
                             front=f$obs[i],
                             rear=m[unitPick,]$obs[mpick],
                             crew=paste(f$obs[i], m[unitPick,][grpPick,]$obs[mpick], sep=""), stringsAsFactors = FALSE
          )
          matches=rbind(matches, newline)
          f$matched[i] = 1
          r$matched[which(pick == TRUE)][which(unitPick == TRUE)][which(grpPick == TRUE)][mpick] = 1
        } else{ #if sum(grpPick) == 0
          #small "open" flocks
          threshold <-  m$grp[unitPick] %in% c(3:5) & f$grp[i] %in% c(3:5)
          diff <- m$grp[unitPick][threshold] - f$grp[i]
          if( sum(threshold) > 0 & f$unit[i] == "open" ){
            mpick <- which( diff == min(diff))[1]
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= (f$grp[i] + m$grp[unitPick][mpick])/2, 
                               unit=f$unit[i],
                               front=f$obs[i],
                               rear=m[unitPick,]$obs[mpick],
                               crew=paste(f$obs[i], m[unitPick,][threshold,]$obs[mpick], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(unitPick == TRUE)][which(threshold == TRUE)][mpick] = 1
          } 
          #large flocks
          threshold <-  m$grp[unitPick] > 5 & f$grp[i] > 5
          diff <- m$grp[unitPick][threshold] - f$grp[i]
          if( sum(threshold) > 0 ){
            mpick <- which( diff == min(diff))[1]
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= (f$grp[i] + m$grp[unitPick][mpick])/2, 
                               unit=f$unit[i],
                               front=f$obs[i],
                               rear=m[unitPick,]$obs[mpick],
                               crew=paste(f$obs[i], m[unitPick,][threshold,]$obs[mpick], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(unitPick == TRUE)][which(threshold == TRUE)][mpick] = 1
          } 
        } #sum(grpPick) > 0
      } #sum(unitPick > 0)  
      if( sum(unitPick) == 0 ){
        #do open/singles
        openPick1 <- (m$unit=="open" & f$unit[i]=="single") | (m$unit=="single" & f$unit[i]=="open")
        openPick2 <- openPick1 & m$grp==f$grp[i] 
        if( sum(openPick2) > 0 ){
          if( !is.na(f$long[1]) ){
            sdist2 = sapply(1:(dim(m[openPick2,])[1]), function(x){
              hav.dist(long1=f[i,"long"], lat1=f[i,"lat"], long2=m[openPick2,][x,"long"], 
                       lat2=m[openPick2,][x,"lat"])}) #in meters
            tdist2 = sapply(1:(dim(m[openPick2,])[1]), function(x){f[i,"ctime"] - m[openPick2,][x, "ctime"]}) #in seconds
            wdist = sqrt( (sdist2/45)^2 + tdist2^2 ) #assume moving at 45 meters per second
          } else{
            wdist = sapply(1:(dim(m[openPick2,])[1]), function(x){f[i,"ctime"] - m[openPick2,][x, "ctime"]}) #in seconds
          }
          mpick = which(wdist==min(wdist))[1]
          newline=data.frame(yr=f$yr[i],
                             tran=f$tran[i],
                             ch="11",
                             sppn=f$sppn[i],
                             grp= (f$grp[i] + m[openPick2,]$grp[mpick])/2,
                             unit=f$unit[i],
                             front=f$obs[i],
                             rear=m[openPick2,]$obs[mpick],
                             crew=paste(f$obs[i], m[openPick2,]$obs[mpick], sep=""), stringsAsFactors = FALSE
          )
          matches=rbind(matches, newline)
          f$matched[i] = 1
          r$matched[which(pick == TRUE)][which(openPick2 == TRUE)][mpick] = 1
        }
        #match front "2 open" with rear "1 pair" or multiple singles
        if(f$grp[i] == 2 & f$unit[i] == "open"){
          pairPick <- m$unit == "pair" & m$grp == 1
          if( sum(pairPick) > 0 ){
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= 1,
                               unit="pair",
                               front=f$obs[i],
                               rear=m[pairPick,]$obs[1],
                               crew=paste(f$obs[i], m[pairPick,]$obs[1], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(pairPick == TRUE)][1] = 1
          } else{
            singlePick <- m$unit == "single" & m$grp == 1
            if( sum(singlePick) >= 2 ){ #make one match of 2 singles
              newline=data.frame(yr=f$yr[i],
                                 tran=f$tran[i],
                                 ch="11",
                                 sppn=f$sppn[i],
                                 grp= 1,
                                 unit="pair",
                                 front=f$obs[i],
                                 rear=m[singlePick,]$obs[1],
                                 crew=paste(f$obs[i], m[singlePick,]$obs[1], sep=""), stringsAsFactors = FALSE
              )
              matches=rbind(matches, newline)
              f$matched[i] = 1
              r$matched[which(pick == TRUE)][which(singlePick == TRUE)][1:2] = 1
            }
          }
        }
        #match rear "2 open" with front "1 pair" or multiple singles
        if(f$grp[i] == 1 & f$unit[i] == "pair"){
          pairPick <- m$unit == "open" & m$grp == 2
          if( sum(pairPick) > 0 ){
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= 1,
                               unit="pair",
                               front=f$obs[i],
                               rear=m[pairPick,]$obs[1],
                               crew=paste(f$obs[i], m[pairPick,]$obs[1], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(pairPick == TRUE)][1] = 1
          } else{
            singlePick <- m$unit == "single" & m$grp == 1
            if( sum(singlePick) >= 2 ){ #make one match of 2 singles
              newline=data.frame(yr=f$yr[i],
                                 tran=f$tran[i],
                                 ch="11",
                                 sppn=f$sppn[i],
                                 grp= 1,
                                 unit="pair",
                                 front=f$obs[i],
                                 rear=m[singlePick,]$obs[1],
                                 crew=paste(f$obs[i], m[singlePick,]$obs[1], sep=""), stringsAsFactors = FALSE
              )
              matches=rbind(matches, newline)
              f$matched[i] = 1
              r$matched[which(pick == TRUE)][which(singlePick == TRUE)][1:2] = 1
            }
          }
        }
        if(f$unit[i] == "flkdrake"){ #match flkdrake to singles and open
          test <- m$unit %in% c("open", "single") & m$grp == f$grp[i]
          if( sum(test) > 0 ){
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= f$grp[i],
                               unit=f$unit[i],
                               front=f$obs[i],
                               rear=m[pairPick,]$obs[1],
                               crew=paste(f$obs[i], m[pairPick,]$obs[1], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(test == TRUE)][1] = 1
          } else{
            singlePick <- m$unit %in% c("open", "single") & m$grp == 1
            if( sum(singlePick) >= f$grp[i] ){ #make a match from mulitple singles
              newline=data.frame(yr=f$yr[i],
                                 tran=f$tran[i],
                                 ch="11",
                                 sppn=f$sppn[i],
                                 grp= f$grp[i],
                                 unit=f$unit[i],
                                 front=f$obs[i],
                                 rear=m[singlePick,]$obs[1],
                                 crew=paste(f$obs[i], m[singlePick,]$obs[1], sep=""), stringsAsFactors = FALSE
              )
              matches=rbind(matches, newline)
              f$matched[i] = 1
              r$matched[which(pick == TRUE)][which(singlePick == TRUE)][1:f$grp[i]] = 1
            }
          }
        } #end match 2 flkdrake to singles
        #what about rear sees flkdrake and front sees single?  How to match that. 
        #match for grp sizes >= 3 for opens and pairs
        #grp size does not need exact match, however must be within "small" (3-5) and large (>5) 
        #flock thresholds
        #small flocks 3-5
        if(f$grp[i] %in% c(3:5) & f$unit[i] == "open"){ # start front observer sees "open"
            #rear sees 2 pair
            test1 <- m$unit == "pair" & m$grp == 2
            if( sum(test1) > 0 ){
              newline=data.frame(yr=f$yr[i],
                                 tran=f$tran[i],
                                 ch="11",
                                 sppn=f$sppn[i],
                                 grp= (f$grp[i] + 2*m[test1,]$grp[1])/2,
                                 unit="open",
                                 front=f$obs[i],
                                 rear=m[test1,]$obs[1],
                                 crew=paste(f$obs[i], m[test1,]$obs[1], sep=""), stringsAsFactors = FALSE
              )
              matches=rbind(matches, newline)
              f$matched[i] = 1
              r$matched[which(pick == TRUE)][which(test1 == TRUE)][1] = 1
            }
            test1 <- m$unit == "flkdrake" & m$grp %in% c(3:5)
            if( sum(test1) > 0 ){
              newline=data.frame(yr=f$yr[i],
                                 tran=f$tran[i],
                                 ch="11",
                                 sppn=f$sppn[i],
                                 grp= (f$grp[i] + 2*m[test1,]$grp[1])/2,
                                 unit="open",
                                 front=f$obs[i],
                                 rear=m[test1,]$obs[1],
                                 crew=paste(f$obs[i], m[test1,]$obs[1], sep=""), stringsAsFactors = FALSE
              )
              matches=rbind(matches, newline)
              f$matched[i] = 1
              r$matched[which(pick == TRUE)][which(test1 == TRUE)][1] = 1
            }
        } # end front observer sees "open"
        if(f$grp[i] == 2 & f$unit[i] == "pair"){ # start front observer sees "2 pair"
          #rear sees 3-5 open
          test1 <- m$unit == "open" & m$grp %in% c(3:5)
          if( sum(test1) > 0 ){
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= (2*f$grp[i] + m[test1,]$grp[1])/2,
                               unit="open",
                               front=f$obs[i],
                               rear=m[test1,]$obs[1],
                               crew=paste(f$obs[i], m[test1,]$obs[1], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(test1 == TRUE)][1] = 1
          }
        } # end front observer sees "2 pair"
        #Large flocks > 5
        if(f$grp[i] > 5 & f$unit[i] == "open") { #front observer sees "open"
          test1 <- m$unit == "pair" & m$grp > 2
          if( sum(test1) > 0 ){
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= (f$grp[i] + 2*m[test1,]$grp[1])/2,
                               unit="open",
                               front=f$obs[i],
                               rear=m[test1,]$obs[1],
                               crew=paste(f$obs[i], m[test1,]$obs[1], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(test1 == TRUE)][1] = 1
          }
        } #front observer sees "open"
        if(f$grp[i] > 2 & f$unit[i] == "pair") { #front observer sees "pair"
          test1 <- m$unit == "open" & m$grp > 5
          if( sum(test1) > 0 ){
            newline=data.frame(yr=f$yr[i],
                               tran=f$tran[i],
                               ch="11",
                               sppn=f$sppn[i],
                               grp= (2*f$grp[i] + m[test1,]$grp[1])/2,
                               unit="open",
                               front=f$obs[i],
                               rear=m[test1,]$obs[1],
                               crew=paste(f$obs[i], m[test1,]$obs[1], sep=""), stringsAsFactors = FALSE
            )
            matches=rbind(matches, newline)
            f$matched[i] = 1
            r$matched[which(pick == TRUE)][which(test1 == TRUE)][1] = 1
          }
        } #front observer sees "pair"
        #do not matched
        if( f$matched[i] != 1){
          newline=data.frame(yr=f$yr[i],
                             tran=f$tran[i],
                             ch="10",
                             sppn=f$sppn[i],
                             grp=f$grp[i],
                             unit=f$unit[i],
                             front=f$obs[i],
                             rear=m$obs[1],
                             crew=paste(f$obs[i], m$obs[1], sep=""), stringsAsFactors = FALSE
          )
          matches=rbind(matches, newline)
          f$matched[i] = 1
        }
      }
    } #length(m) >0
    
  } #i loop
  

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