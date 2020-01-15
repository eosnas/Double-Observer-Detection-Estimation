#for 2015 and 2016 double observer data sets
#Read in transcribed data, format to make usable, and then pass to match function
#code written by Erik Osnas, 7 July 2018

library(readxl)
library(RCurl)

#read in QA/QC'ed data files
front <- as.data.frame(read_excel(
  path="K:/ACP SPEI Detection Survey/Clean Sheet Final Data/acp2 double obs front seat.xls",
  sheet=1, col_names = TRUE, 
  col_types = c("numeric", "text", "skip", "numeric", "numeric", "numeric", "numeric", 
                "text", "text", "text", "numeric", "text", rep("text", 7))))

rear <- as.data.frame(read_excel(
  path="K:/ACP SPEI Detection Survey/Clean Sheet Final Data/acp2 double obs rear seat.xls",
  sheet=1, col_names = TRUE, 
  col_types = c("numeric", "text", "skip", "numeric", "numeric", "numeric", "numeric", 
                "text", "text", "text", "numeric", "text", rep("text", 7))))
front <- front[is.na(front$drop),]
rear <- rear[is.na(rear$drop),]
str(front)
str(rear)

#remove unwanted species, hens, begin and end points:
#Note that the file above are not "Green light" files
spplist <- c("GWFG",  "PALO",  "KIEI",  "GLGU",  "NOPI",  "LTDU",  "SAGU",  "TUSW",  "ARTE",  "JAEG", 
             "YBLO",  "CAGO",  "RTLO",  "SCAU",  "SPEI",  "SNGO", "UNSB",  "BLBR")
front <- front[front$sppn %in% spplist,]
rear <- rear[rear$sppn %in% spplist,]
#convert "num" to numeric
front$num <- as.numeric(front$num)
rear$num <- as.numeric(rear$num)
str(front)
str(rear)
#change names for num and grp
names(front)[12:13] <- names(rear)[12:13] <- c("grp", "unit")
front$obsvr <- toupper(front$obsvr)
rear$obsvr <- toupper(rear$obsvr)

#Add NA log and lat columns
front$long <- front$lat <- rep(NA, dim(front)[1])
rear$long <- rear$lat <- rep(NA, dim(rear)[1])

#inspect data
unique(front$trans)
unique(rear$trans)
prop.table(table(front$unit, front$obsvr))
prop.table(table(rear$unit, rear$obsvr))

source("findNear.R")
#test DubMatch function
# #make data subset
#front <- front[front$sppn == "SPEI",]
#rear <- rear[rear$sppn == "SPEI",]
# for(i in 1:dim(front)[1]){
#   s <- findNear(st=as.numeric(front[i,c("long", "lat", "adjctime", "day", "mo", "year")]),
#                 y=rear[,c("long", "lat", "adjctime", "day", "mo", "year")],
#                 sdist=250, tdist=5, tlag=0)
#   print(front[i,"obsvr"])
#   print(front[i,c("seat", "adjctime", "sppn", "grp", "unit")])
#   print(rear[s,c("obsvr", "seat", "adjctime", "sppn", "grp", "unit")])
# }

source("DubMatch2.R")
dat <- rbind(front, rear)
names(dat) <- c("finalorder", "drop", "seq", "yr", "mo", "da", "se", "obs", "tran", "ctime",
                "sppn", "grp", "unit", "distance", "easy", "behav", "other", "notes", "lat", "long")
datr <- dat[dat$obs == "WWL" & dat$se == "lr",]
trans <- unique(datr$tran)
datf <- dat[dat$obs == "HMW" & dat$tran %in% trans & dat$yr == 2016,]  #Bill sat in back only in 2016
match1 <- DubMatch(front.data = datf, rear.data=datr, front.obs = "HMW", rear.obs = "WWL")
datr <- dat[dat$obs == "TKZ" & dat$se == "lr",]
trans <- unique(datr$tran)
datf <- dat[dat$obs == "HMW" & dat$tran %in% trans,]
match2 <- DubMatch(front.data = datf, rear.data=datr, front.obs = "HMW", rear.obs = "TKZ")
datr <- dat[dat$obs == "TKZ" & dat$se == "rr",]
trans <- unique(datr$tran)
datf <- dat[dat$obs == "WWL" & dat$tran %in% trans,]
match3 <- DubMatch(front.data = datf, rear.data=datr, front.obs = "WWL", rear.obs = "TKZ")
match <- rbind(match1, match2, match3)
match[match$sppn == "SPEI",]

#table all matches by spp
table(match$ch, match$sppn)
table(match$sppn)
match$grpcode <- ifelse(match$grp == 1 & match$unit %in% c("open", "single"), 1, 
                        ifelse(match$grp == 1 & match$unit == "pair", 2, 
                               ifelse(match$grp == 2 & match$unit %in% c("open", "flkdrake", "single"), 2, 
                                      ifelse(match$grp >= 3 & match$grp <= 5 & match$unit %in% c("open", "flkdrake", "single"), 3, 
                                             ifelse(match$grp == 2 & match$unit == "pair", 3, 
                                                    ifelse(match$grp > 5 & match$unit %in% c("open", "flkdrake", "single"), 4, 
                                                           ifelse(match$grp >= 3 & match$unit == "pair", 4, "other")
                                                           )
                                                    )
                                             )
                                      )
                               )
                        )


