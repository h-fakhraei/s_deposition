
# Code for the kriging analysis and modelling of atmospheric S deposition 

# Description --------------------------------------------------------- 
# This repository includes the kriging model and input files used for E-L. S. Hinckley,  
# J. T. Crawford1, H. Fakhraei and C. T. Driscoll for A shift in sulfur-cycle manipulation from
# atmospheric emissions to agricultural additions, Nature Geoscience 2020. This is a kriging model to 
# interpolate atmospheric deposition through no measured locations throughout the US using point observation
# from two monitoring networks (i.e., NADP and CASTNET)
# 
# You are welcome to use or adapt this code providing you cite
# E-L.S. Hinckley, J. T. Crawford1, H. Fakhraei and C.T. Driscoll for A shift in sulfur-cycle manipulation from
# atmospheric emissions to agricultural additions, Nature Geoscience 2020.
# 
# To run this code a set of input data is required to be downloaded from the public domain resources.
# Here, we included these datasets until year 2017 but the user will be able to download the data and 
# perform kriging as the datasets are getting updated.

# Inputs --------------------------------------------------------------- 
# NADP wet deposition data (http://nadp.slh.wisc.edu/ntn/)
# location of wet deposition sites (http://nadp.slh.wisc.edu/data/sites/NTN/?net=NTN)
# CASTNET dry deposition data (https://java.epa.gov/castnet/clearsession.do)
# location of dry deposition sites (https://www.epa.gov/castnet/castnet-site-locations)
# PRISM precipitation quantity (ftp://prism.oregonstate.edu/monthly/ppt/) 

# Outputs --------------------------------------------------------------- 
# For each running year the code generates three GeoTiff maps including wet, dry and total S deposition for the entire US.
# as example output maps for years 1989 and 2017 are included here.

# Required package ------------------------------------------------------ 
library(geoR)
library(fields)
library(maps)
library(rgdal)
library(raster) 

# load wet deposition data and NADP site locations----------------------- 
# download wet deposition data from NADP website: http://nadp.slh.wisc.edu/ntn/
# note units are in mg/l and e.g., sulfate is in SO4 mg/l not in S mg/l
nadp_data<-read.csv("NTN-All-cy.csv")
# or this can be loaded directly from the web by:
#nadp_data<-read.csv("http://nadp.slh.wisc.edu/datalib/ntn/cy/NTN-All-cy.csv")

nadp_data[nadp_data==-9] <- NA
names(nadp_data)[names(nadp_data)=="siteID"] <- "SiteID"

# only include years that have been sampled for at least 330 days
nadp_data<-nadp_data[nadp_data$daysSample>330,]

# download location of wet deposition sites from: http://nadp.slh.wisc.edu/data/sites/NTN/?net=NTN
nadp_loc<-read.csv("NTNsites.csv",na.strings = "")
# or this can be loaded directly from this address:
#nadp_loc<-read.csv("http://nadp.slh.wisc.edu/data/sites/CSV/?net=NTN",na.strings = "")
names(nadp_loc)[2] <- "SiteID"

# remove Hawaii and Alaska (just focus on the continental states)
nadp_loc <- nadp_loc[!(nadp_loc$state%in%c("HI","AK")),] 

# combine wet deposition data and locations
nadp <-merge(nadp_loc,nadp_data,by="SiteID")


# load dry deposition data and CASTNET site locations---------------------- 
# download annual dry deposition data from CASTNET website: https://java.epa.gov/castnet/clearsession.do
# following csv file downloaded from CASTNET website on 11/4/2018 
# note units are in kg/ha and e.g., sulfate is in kg SO4/ha not in Kg S/ha
# note recently the format and unit of dry deposition data on CASTNET website have been changed
castnet_data <-read.csv("Dry_Deposition_Annual_all_US_kgha_2017.csv",na.strings="")

# rename castnet data frame
names(castnet_data)[names(castnet_data)=="SITE_ID"] <- "SiteID"
names(castnet_data)[names(castnet_data)=="YEAR"] <- "yr"
names(castnet_data)[names(castnet_data)=="SO2_FLUX"] <- "SO2"
names(castnet_data)[names(castnet_data)=="SO4_FLUX"] <- "SO4"

# remove CASTNET locations without SO2 record
castnet_data<-castnet_data[!is.na(castnet_data$SO2),]

# order castnet data
castnet_data<-castnet_data[order(castnet_data$SiteID,castnet_data$yr),]

# download location of dry deposition sites from CASTNET website:
# https://www.epa.gov/castnet/castnet-site-locations
castnet_loc <-read.csv("castnet_locations_entire_US.csv",na.strings="")
names(castnet_loc)[names(castnet_loc)=="SITE_ID"] <- "SiteID"
castnet_loc <- castnet_loc[order(castnet_loc$SiteID),]
# change column names
names(castnet_loc)[names(castnet_loc)=="Latitude"] <- "latitude"
names(castnet_loc)[names(castnet_loc)=="Longitude"] <- "longitude"
names(castnet_loc)[names(castnet_loc)=="Elevation"] <- "elevation"

# Remove Hawaii and Alaska (just focus on the continental states)
castnet_loc <- castnet_loc[!(castnet_loc$STATE%in%c("HI","AK")),] 

# combine dry deposition data and locations
castnet<-merge(castnet_loc,castnet_data,by="SiteID")

# consider S dry deposition as dry SO2 + dry SO4
castnet$SO4 <- castnet$SO4 + castnet$SO2*96.06/64.06
castnet$SO2 <- NULL


# Ordinary Kriging----------------------------------------------------------
# define boundaries for the kriging (to use continental sites)
max_lat<-49.93750;min_lat<-24.0625;max_long<-(-66.4791667);min_long<-(-125.0208333)

# define the kriging function
do_krig <- function(d,s,model){
  
  # Variogram analysis
  up<-min(abs(max_lat-min_lat),abs(min_long-max_long))
  bins<-seq(0,up,up/(length(d)/5))  #there are the bin edges
  vg <- variog(data=d,coords=s)
  plot(vg,pts.range = c(1,3))
  
  # fit model to vario :
  ml.exp <- variofit(vg)
  ml.mat1.5  <- variofit(vg, kappa=1.5)
  ml.sph  <- variofit(vg, cov.model="sph")
  ml.gau  <- variofit(vg, cov.model = "gau")
  ml.nug  <- variofit(vg, cov.model = "pure.nugget")
  
  lines(ml.exp,col="red")
  lines(ml.mat1.5,col="blue")
  lines(ml.sph,col="green")
  lines(ml.gau,col="black")
  lines(ml.nug,col="gray",lty=2,lwd=2)
  legend("topleft",c("Exponential","Matern","Spherical","Gaussian","Pure Nugget"),
         col=c("red","blue","green","black","gray"),lty=c(1,1,1,1,2),
         lwd=c(1,1,1,1,2))
  
  # Create grid of prediction points based on PRISM grids
  cellsize <- 0.04166667
  sp1<-seq(min_long,max_long,by=cellsize)
  sp2<-seq(min_lat,max_lat,by=cellsize)
  
  sp<-expand.grid(sp1,sp2)
  
  # if there is an issue with one model type (exp) another model being used :
  # first try the requested model by user  
  # then if it does not work try another one in the "models" vector 
  models <- c("exp","mat","sph","gau","nug")
  for (modtype in c(model,models[models!=model])){
    if(modtype=="exp") mod<-ml.exp
    if(modtype=="mat") mod<-ml.mat1.5
    if(modtype=="sph") mod<-ml.sph
    if(modtype=="gau") mod<-ml.gau
    if(modtype=="nug") mod<-ml.nug#Pure Nugget
    pred<-try(krige.conv(data=d,coords=s,locations=sp,krige=krige.control(obj.m = mod)))
    if(class(pred)!="try-error") break
  }
  
  # Plot the predicted values:
  image.plot(sp1,sp2,matrix(pred$predict,length(sp1),length(sp2)),
             zlim=range(d),ylim=c(min_lat,max_lat),xlim=c(min_long,max_long),
             las=1,ylab="Latitude",xlab="Longitude")
  #contour(pred, add=TRUE, nlev=21)
  
  rotate <- function(x) t(apply(x, 2, rev))
  xy <- matrix(pred$predict,length(sp1),length(sp2))
  xy<- rotate(rotate(rotate(xy)))
  return(xy)
}

# define initial year for starting the kriging
initial_year<-1989 # CASTNET data started from 1987 but enough data for kriging since 1989
# define the variogram model 
model="exp"# variogram model other options are:  c("mat","sph","gau","nug")

# run kriging for each individual year starting 1989
for (year in initial_year:2017){
  
  # extract and rearrange the required data from NADP file:
  # conc. of wet deposition:
  nadp2<-data.frame(SiteID=nadp$SiteID,LATITUDE=nadp$latitude,
                    LONGITUDE=nadp$longitude,ELEVATION=nadp$elevation,
                    Year=nadp$yr,param=nadp[,"SO4"])
  nadp2<-subset(nadp2,nadp2$LATITUDE>min_lat & nadp2$LATITUDE<max_lat
                &nadp2$LONGITUDE>min_long & nadp2$LONGITUDE<max_long)
  nadp2<-subset(nadp2,nadp2$Year==year)
  nadp2<-na.omit(nadp2)
  dat_nadp<-nadp2[,6]#observation
  loc_nadp<-cbind(nadp2[,3],nadp2[,2])#lat,long
  
  # extract and rearrange the required data from CASTNET file:
  # flux of dry deposition:
  castnet2<-data.frame(SiteID=castnet$SiteID,LATITUDE=castnet$latitude,
                       LONGITUDE=castnet$longitude,ELEVATION=castnet$elevation,
                       Year=castnet$yr,param=castnet[,"SO4"])
  castnet2<-subset(castnet2,castnet2$LATITUDE>min_lat & castnet2$LATITUDE<max_lat
                   &castnet2$LONGITUDE>min_long & castnet2$LONGITUDE<max_long)
  castnet2<-subset(castnet2,castnet2$Year==year)
  castnet2<-na.omit(castnet2)
  dat_cast<-castnet2[,6]#observation
  loc_cast<-cbind(castnet2[,3],castnet2[,2])#lat,long
  
  # conduct kriging
  k_nadp<- do_krig (dat_nadp,loc_nadp, model)
  k_cast<- do_krig (dat_cast,loc_cast, model)
  
  # extract precipitation data from already downloaded data from the PRISM model
  # data downloaded from: ftp://prism.oregonstate.edu/monthly/ppt/
  # first unzip the "PRISM" files then load the precipitation data
  # Note because of size limit for uploading the PRISM data files on Github they are split into two folders:
  ifelse(year>2000,prismfolder<-"PRISM_2001_2017",prismfolder<-"PRISM_1989_2000")
  setwd(paste0("./",prismfolder,"/PRISM_ppt_stable_4kmM3_",year,"_all_bil"))
  ss<-readGDAL(paste0("PRISM_ppt_stable_4kmM3_",year,"_bil.bil"))
  prism<-t(matrix(ss@data[[1]],1405,621))
  
  # calculate wet dep= precipitation* conc.
  wdep <- prism * k_nadp /100 * 32.06/96.06# prism: mm, conc: mgl -->kg SO4/ha-->Kg S/ha
  # note above code caused NA for ooutside US land (since prism is NA for those cells)
  
  # to make castnet cells outside us land NA do this:
  ddep <- k_cast *prism/prism
  ddep <- ddep * 32.06/96.06# e.g., SO4 to S kg/ha 
  tdep <- wdep + ddep
  
  # Turn the matrix into a raster
  rast_tdep <- raster(tdep) 
  rast_ddep <- raster(ddep) 
  rast_wdep <- raster(wdep)
  # define the extents
  extent(rast_tdep) <- c(min_long,max_long,min_lat,max_lat)
  extent(rast_ddep) <- c(min_long,max_long,min_lat,max_lat)
  extent(rast_wdep) <- c(min_long,max_long,min_lat,max_lat)
  # ... and assign a projection (the same as PRISM)
  crs(rast_tdep) <- crs(ss) 
  crs(rast_ddep) <- crs(ss)
  crs(rast_wdep) <- crs(ss)
  
  # generate Kriging maps in GeoTiff format
  # unit: Kg S/ha (not Kg SO4/ha):
  setwd('..')
  setwd('..')
  setwd("./output_maps")
  writeRaster(rast_tdep, paste0("./TDEP_",year), format = "GTiff",overwrite=TRUE,prj=T)
  writeRaster(rast_ddep, paste0("./DDEP_",year), format = "GTiff",overwrite=TRUE,prj=T)
  writeRaster(rast_wdep, paste0("./WDEP_",year), format = "GTiff",overwrite=TRUE,prj=T)
  setwd('..')
}      


