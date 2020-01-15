#' Find rows of a dataframe that are close in time and/or space (matches)
#' @param st A vector that gives spatial and temporal location (Long, Lat, seconds, day, month, year).
#' @param y A dataframe with columns to match st.  
#' @param sdist Spatial distance to determine match.
#' @param tdist Temporal distance to determine match.
#' @param tlag Temporal lag applied to y.
#' @return A logical vector that indicates matches in y given other parameters. 
#'
#' @examples
findNear <- function(st=NA, y=NA, sdist=250, tdist=5, tlag=0){
  if( !is.na(st[1]) ){
    space <- sapply(1:(dim(y)[1]), function(x){hav.dist(long1=st[1], lat1=st[2], long2=y[x,1], lat2=y[x,2])})
    time <- abs(st[3]-y[,3])
    match <- space < sdist & st[4]==y[,4] & st[5]==y[,5] & st[6]==y[,6] & time < (tdist + tlag)
  } else{
    time <- abs(st[3]-y[,3])
    match <- st[4]==y[,4] & st[5]==y[,5] & st[6]==y[,6] & time < (tdist + tlag)
  }
  return(match)
}

#' Calculate Haversine distance based on Longitude and Latitude of two points
#'
#' @param long1 
#' @param lat1 
#' @param long2 
#' @param lat2 
#'
#' @return A number in meters.  
#' @export
#'
#' @examples
hav.dist <- function(long1, lat1, long2, lat2) {
  long1=deg2rad(long1)
  lat1=deg2rad(lat1)
  long2=deg2rad(long2)
  lat2=deg2rad(lat2)
  R <- 6371
  diff.long <- (long2 - long1)
  diff.lat <- (lat2 - lat1)
  a <- sin(diff.lat/2)^2 + cos(lat1) * cos(lat2) * sin(diff.long/2)^2
  b <- 2 * asin(pmin(1, sqrt(a))) 
  d = R * b
  return(d)
}

#' Title
#'
#' @param deg 
#'
#' @return
#' @export
#'
#' @examples
deg2rad <- function(deg){ return(deg*pi/180)
}