
# Update R script for 4 areas (QH, 4/24/2026)

#IMPORT DATA
#setwd('C:/users/matthew.lauretta/desktop/2026/bluefin/etag_data_summary/')
#data=read.csv('BFT_geolocations_2026_03_31.csv')
data=read.csv('data/Etag/BFT_geolocations_2026_03_31.csv')
head(data)

#DEFINE DATA FIELDS REFERENCES FOR THE DATASET TO MATCH THE CODE BELOW
lat = data$lat+0.000001
lon = data$lon+0.000001

data$Reference_ID = data$tag
data$Date = data$date

#REQUIRED R PACKAGES
	library(sp)
	library(maps)
	library(mapdata)

windows()
	maps::map('worldHires',col=c('gray'),fill=T,xlim=c(-100,45),ylim=c(-50,80))
	axis(1,at=seq(-100,45,5))
	axis(2,at=seq(-50,80,5))
	mtext('Longitude',1,line=3)
	mtext('Latitude',2,line=3)
points(data$lon,data$lat,pch=".",cex=2,col=2)


## FILL OUT THE MISSING SIZES FROM THE WEIGHT
tags <- unique(data$tag[which(is.na(data$size_at_tagging)==TRUE)])
for (i in 1:length(tags)) {
     ind <- which(data$tag==tags[i])
     if (is.na(data$wt_at_tagging[ind[1]])==FALSE) {
         data$size_at_tagging[ind] <-
trunc((as.numeric(data$wt_at_tagging[ind[1]])/0.0000315551)^(1/2.8984539))
     }
}
tags2 <- unique(data$tag[which(is.na(data$size_at_tagging)==TRUE)])

#DEFINED STOCK AREA X (LON) AND Y (LAT) BOUNDARIES
	BFT1=list(x=c(-80,-88,-95,-100,-100,-85,-80), y=c(20,20,16.5,20,35,35,25))
	BFT2=list(x=c(-82.5,-75,-75,-65,-65,-55,-55,-70,-95,-88,-80,-80,-82.5),
		y=c(30,30,25,25,20,20,0,0,16.5,20,20,25,30))
	BFT3=list(x=c(-70,-70,-60,-55,-55), y=c(45,55,55,50,45))
	BFT4=list(x=c(-70,-55,-55,-65,-65,-75,-75,-82.5,-85,-70,-55,-55,-60,-70,-80,-100,-100,-45,-45,-30,-30,-25,-25,-70),
		y=c(0,0,20,20,25,25,30,30,35,45,45,50,55,55,50,60,80,80,10,10,5,5,-50,-50))
	BFT5=list(x=c(-30,-45,-45,-30), y=c(40,40,80,80))
	BFT6=list(x=c(-30,-45,-45,-30), y=c(10,10,40,40))
	BFT7=list(x=c(-30,45,45,15,15,-15,-15,-30,-30), y=c(80,80,50,50,60,60,50,50,80))
	BFT8=list(x=c(-30,-30,-15,-15,15,15,5,-5), y=c(40,50,50,60,60,50,50,40))
	BFT9=list(x=c(-30,-30,-5,-5,20,20,-25,-25,-30), y=c(10,40,40,30,30,-50,-50,5,5))
	BFT10=list(x=c(-5,-5,5,23,23), y=c(30,40,50,50,30))
	BFT11=list(x=c(23,45,45,23) ,y=c(50,50,30,30))

#STOCK AREAS PLOTTED ON MAP
	maps::map('worldHires',col=c('gray'),fill=T,xlim=c(-100,45),ylim=c(-50,80))
	axis(1,at=seq(-100,45,5))
	axis(2,at=seq(-50,80,5))
	mtext('Longitude',1,line=3)
	mtext('Latitude',2,line=3)
	polygon(BFT1,border=1,lwd=2)
	text("GOM",x=-90,y=25,font=2,col=2)
	polygon(BFT2,border=1,lwd=2)
	text("WATL",x=-70,y=15,font=2,col=2)
	polygon(BFT3,border=1,lwd=2)
	text("WATL",x=-63,y=49,font=2,col=2)
	polygon(BFT4,border=1,lwd=2)
	text("WATL",x=-60,y=30,font=2,col=2)
	polygon(BFT5,border=1,lwd=2)
	text("EATL",x=-37.5,y=60,font=2,col=2)
	polygon(BFT6,border=1,lwd=2)
	text("EATL",x=-37.5,y=25,font=2,col=2)
	polygon(BFT7,border=1,lwd=2)
	text("EATL",x=-5,y=70,font=2,col=2)
	polygon(BFT8,border=1,lwd=2)
	text("EATL",x=-10,y=45,font=2,col=2)
	polygon(BFT9,border=1,lwd=2)
	text("EATL",x=-5,y=0,font=2,col=2)
	polygon(BFT10,border=1,lwd=2)
	text("MED",x=10,y=40,font=2,col=2)
	polygon(BFT11,border=1,lwd=2)
	text("MED",x=35,y=40,font=2,col=2)

#STOCK AREA ASSIGNMENT BASED ON BFT LOCATION IN DECIMAL DEGRESS LAT AND LON
	BFT_area=c("GOM", "WATL", "WATL", "WATL", "EATL", "EATL", "EATL", "EATL", "EATL", "MED", "MED")
	data$STOCK_AREA=as.character(sapply(1:length(data[,1]),function(i)BFT_area[which(c(
		point.in.polygon(lon[i],lat[i],BFT1$x,BFT1$y),
		point.in.polygon(lon[i],lat[i],BFT2$x,BFT2$y),
		point.in.polygon(lon[i],lat[i],BFT3$x,BFT3$y),
		point.in.polygon(lon[i],lat[i],BFT4$x,BFT4$y),
		point.in.polygon(lon[i],lat[i],BFT5$x,BFT5$y),
		point.in.polygon(lon[i],lat[i],BFT6$x,BFT6$y),
		point.in.polygon(lon[i],lat[i],BFT7$x,BFT7$y),
		point.in.polygon(lon[i],lat[i],BFT8$x,BFT8$y),
		point.in.polygon(lon[i],lat[i],BFT9$x,BFT9$y),
		point.in.polygon(lon[i],lat[i],BFT10$x,BFT10$y),
		point.in.polygon(lon[i],lat[i],BFT11$x,BFT11$y))==1)]))
write.csv(data,'data/Etag/BFT_geolocations_areas.csv',row.names=FALSE)

#DATA AGGREGATION: DAYS PER STOCK AREA BY TAG_ID AND REGION TRANSITION
	data$REGION_ENTRY=1
	for(i in 2:length(data[,1]))
		{
		data$REGION_ENTRY[i]=ifelse(data$Reference_ID[i]==data$Reference_ID[i-1]&data$STOCK_AREA[i]==data$STOCK_AREA[i-1],data$REGION_ENTRY[i-1],data$REGION_ENTRY[i-1]+1)
		}
	summary=aggregate(data$Date,by=list(data$group,data$Reference_ID,data$size_at_tagging,data$STOCK_AREA,data$REGION_ENTRY),length)
	colnames(summary)=c('Group','Tag_ID','Size_at_release','Stock_Area','Entry','Days')
	summary$Start_Date=sapply(1:length(summary[,1]),function(i)data$Date[data$Reference_ID==summary$Tag_ID[i]&
		data$STOCK_AREA==summary$Stock_Area[i]&data$REGION_ENTRY==summary$Entry[i]][1])
	summary$End_Date=sapply(1:length(summary[,1]),function(i)rev(data$Date[data$Reference_ID==summary$Tag_ID[i]&
		data$STOCK_AREA==summary$Stock_Area[i]&data$REGION_ENTRY==summary$Entry[i]])[1])
	#summary
write.csv(summary,'data/Etag/BFT_etags_processed_forM3_20260424.csv',row.names=FALSE)




