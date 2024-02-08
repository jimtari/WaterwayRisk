### This file prepares the threat layers for Zonation. Exactly what is needed will
# depend on the threat but the output layer will be a raster of wetlands with 
# a scaled threat value between 0 (no threat) and 1 (maximum threat) 

#Created 14 Nov 2023 - adapted to WIT data based on Prepare_threats_4Zonation.R
# Last updated 21 Nov 2023 - removed CIs, made sure new results sent to 
# C:/Data/AAAMoore/R/Wetland_Indices/PilotV4
#append _WIT to the _raster threat outputs

# CIs were CIs of bootstrapped null model not the CIs of the stress index.  
# UMelb team are going to look into that later but in the meantime the stress index is ok. 

# 15 Nov 2023 - WILL NEED TO RERUN WHEN STRESS DATA UPDATED (CIs ARE WRONG) 
#ALSO CHECK THAT writes to correct place etc.  I've decided to keep in 
# original R project and call WIT stress index PilotV4
# 19 Dec 2023 - Reran with the correct confidence intervals - overworte stuff in PilotV4

library(tidyverse)
library(sf)
library(stars)
library(terra)
library(mapview)

#########################################################################################
## Import stress index and infer missing values - only needs to be done once - output stored as rds
########################################################################################

# stress index

# 1. Characteristic: mean inundation magnitude
# 2. Infer missing wetlands using nearest neighbour of same wetland type
#   i. read in geopackage
#   ii. infer stress based on same value as nearest neighbour of same wetland type - infer mean, LCI and UCI at same time
# 3. Just negative stress values used as this corresponds to associations (but does not address coastal inundation)
#    i. set everything greater than 0 as 0
#    ii. take absolute values to change to positive
# 4. Stress function: Linear for now
#    i. standardise to between zero and one (divide by max value)
#    ii. consider rounding so the numbers aren't so ridiculous.
# 5. rasterise threat layer to standardised number to fit wetland raster and output
# 6. repeat 3-5 using LCI and UCI
# 7. repeat 3-6 without inferring missing data

# #1. read in geopackage
# wetland_stress <- sf::st_read("D:/Wetland Indices Pilot/UMelb Hydrological stress Do NOT CIRCULATE/ForJos_v2_v3_combined/Wetland_stress.gpkg")  %>%
#   select(UID, RecID, WETLAND_NO, HECTARES, IWG_num, IWG_ID, Int_grp, Int_grp_crit,XTNDWT,
#          XWT_COL,SIMPWT,SWT_COL,ORIGIN, TENURE, SIMPTNR, CMA_ID, CMA, CMA_ABBR, Mag,
#          Mag_dt, Mag_Max, Mag_Maxdt, Mag_Min, Mag_Mindt, MagStress, MagSig,
#          MagStressLCI, MagStressUCI) %>%
#         mutate(MagStressInfer = MagStress, MagStressInferLCI = MagStressLCI, MagStressInferUCI = MagStressUCI, .after = MagStressUCI)
# 
# wetland_stress$MagSigClass <- "Not significant"
# wetland_stress$MagSigClass[wetland_stress$MagSig <= 0.05] <- "Significant"
# wetland_stress <- wetland_stress %>%
#   relocate("MagSigClass", .after="MagSig")
# 
# 
# names(wetland_stress)
# infer <- wetland_stress %>%
#   filter(is.na(MagStressInfer))
# 
# Ninfer <- nrow(infer)
# 
# sf::sf_use_s2(TRUE)  #nearest neighbour works better with this set to TRUe I think
# closest <- 0
# 
# # This loop takes about half an hour for ~3000 wetlands
# for(i in 1:Ninfer){
#   wetland1<- subset(wetland_stress, !(wetland_stress$RecID %in% infer$RecID) & wetland_stress$XTNDWT == infer$XTNDWT[i]) #choose wetlands with same type wetland and no missing data
#   closest <- st_nearest_feature(infer[i,], wetland1, check_crs = T, longlat = T)
#   wetland_stress$MagStressInfer[wetland_stress$RecID == infer$RecID[i]] <- wetland1$MagStressInfer[closest]
#   wetland_stress$MagStressInferLCI[wetland_stress$RecID == infer$RecID[i]] <- wetland1$MagStressInferLCI[closest]
#   wetland_stress$MagStressInferUCI[wetland_stress$RecID == infer$RecID[i]] <- wetland1$MagStressInferUCI[closest]
#   print(i)
# }
# 
# sf::sf_use_s2(FALSE)
# 
# # # This is a test to see if working
# # infer$RecID[254]
# # wetland_stress$MagStressInfer[which(wetland_stress$RecID==6754)]
# #
# # wetland_stress$RecID[1509]
# # wetland_stress$MagStressInfer[which(wetland_stress$RecID==24807)]
# 
# saveRDS(wetland_stress, paste0(getwd(), "/PilotV4/wetland_stress.rds"))

##################################################################################
## Process the stress layer, convert to raster and save as tif
##################################################################################         

# 3. Just negative stress values used as this corresponds to associations (but does not address coastal inundation)
#    i. set everything greater than 0 as 0
#    ii. take absolute values to change to positive
# 4. Stress function: Linear for now  
#    i. standardise to between zero and one (divide by max value)
#    ii. consider rounding so the numbers aren't so ridiculous.I rounded to six decimal places
# 5. rasterise threat layer to standardised number to fit wetland raster and output
# 6. repeat 3-5 using LCI and UCI  
##  Not doing this for now
# 7. repeat 3-6 without inferring missing data


 
setwd("O:/SAN_Projects/WaterwaysValuesThreats/Indices_code/Wetland_Indices")

wetland_stress <- readRDS(paste0(getwd(), "/PilotV4/wetland_stress.rds"))

threat <-  "MagStressInfer"  

threatname <-  "HydroStress_dry"  

workingsf <- wetland_stress %>% 
  select(UID, RecID, "threat" = all_of(threat))

#testing section
test1 = nrow(workingsf[workingsf$threat>0 & !is.na(workingsf$threat), ])
testna1 = nrow(workingsf[is.na(workingsf$threat), ])

#NOT rescaling - just raserterize the threat values....
#calc section
#workingsf = workingsf %>%
#  rowwise() %>%
#    mutate(threat = abs(min(threat, 0)))
 
#testing section
max(workingsf$threat, na.rm=T) 
min(workingsf$threat, na.rm=T) 

test2 = nrow(workingsf[workingsf$threat==0 & !is.na(workingsf$threat), ])
testna2=nrow(workingsf[is.na(workingsf$threat), ])
plot(density(workingsf$threat, na.rm=T))

#calc section
maxthreat <- max(workingsf$threat, na.rm=T)

#NOT rescaling - just raserterize the threat values....
#workingsf = workingsf %>%
#  mutate(threat = round(threat/maxthreat,6))   #I am rounding to 6 decimal places here (i think)

#testing section testX should all be same number as should testmaX
max(workingsf$threat, na.rm=T) 
min(workingsf$threat, na.rm=T)   

test3 = nrow(workingsf[workingsf$threat==0 & !is.na(workingsf$threat), ])
testna3 = nrow(workingsf[is.na(workingsf$threat), ])
plot(density(workingsf$threat, na.rm=T))

#rasterise

in_directory = "O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data"
out_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data"

#read in the wetlands tif as a template
raster_template <- read_stars(file.path(in_directory, "wetlands_WIT_raster.tif"))
raster_template[[1]] <- NA  #need to do this as otherwise end up with an imprint of template raster values



wetlandRas<-rast(file.path(in_directory, "wetlands_WIT_raster.tif"))

layer2rasterise <- workingsf %>%
  select(threat) 

# min(layer2rasterise$threat, na.rm=T)
# max(layer2rasterise$threat, na.rm=T)

rast_template=rast(raster_template)

wetland_threat<- terra::rasterize(st_transform(layer2rasterise,crs(rast_template)), rast_template,field='threat',touches=TRUE,fun='min')

# mapview(layer2rasterise)+
#   mapview(threat_stars , downsample=20)

writeRaster(wetland_threat,file.path(out_directory, paste0(threatname, "_raster_wetland.tif")))

##isc reaches..

isc=st_read("O:/SAN_Projects/WaterwaysValuesThreats/Asset_layers/ISC2010_FSR2020_stats.shp")

iscbuff<-st_buffer(isc,dist=75)

iscbuff$ID=as.numeric(gsub("_","",iscbuff$UNIQUEID))

#iscbuff$threat<-(iscbuff$Q90s+iscbuff$PZDs)/2  #wrong

iscbuff$threat<-iscbuff$Q10

layer2rasterise <- iscbuff %>%  select(threat) 


reach_threat <- terra::rasterize(st_transform(layer2rasterise,crs(rast_template)), rast_template,field='threat',touches=TRUE,fun='max')

writeRaster(reach_threat,file.path(out_directory, paste0(threatname, "_raster_reach.tif")),overwrite=T)


### estuary

ests=st_read("O:/SAN_Projects/WaterwaysValuesThreats/Asset_layers/Estuaries_10Oct2023.gpkg")

estuary_threat <- terra::rasterize(st_transform(ests,crs(rast_template)), rast_template,touches=TRUE)*0

writeRaster(estuary_threat,file.path(out_directory, paste0(threatname, "_raster_estuary.tif")))

#end of section that standardises the stress layer

#CI section below#
##################################################################################################
## repeat for LCI and UCI
####  LCI
threat <-  "MagStressInferLCI" #"MagStressInferUCI"

workingsf <- wetland_stress %>% 
  select(UID, RecID, "threat" = all_of(threat))

#testing section
max(workingsf$threat, na.rm=T) 
min(workingsf$threat, na.rm=T) 

#calc section
workingsf = workingsf %>%
  rowwise() %>%
  mutate(threat = abs(min(threat, 0)))

#testing section
max(workingsf$threat, na.rm=T) 
min(workingsf$threat, na.rm=T) 

#calc section
maxthreat <- max(workingsf$threat, na.rm=T)
workingsf = workingsf %>%
  mutate(threat = round(threat/maxthreat,6))   #I am rounding to 6 decimal places here (i think)

#rasterise
in_directory = "D:/Wetland Indices Pilot/Data/"
out_directory <-"D:/Wetland Indices Pilot/Data/"

#read in the wetlands tif as a template
raster_template <- read_stars(paste0(in_directory, "wetlands_WIT_raster.tif"))
raster_template[[1]] <- NA  #need to do this as otherwise end up with an imprint of template raster values

layer2rasterise <- workingsf %>%
  select(threat) 

threat_stars <- st_rasterize(layer2rasterise, template=raster_template, options = "ALL_TOUCHED=TRUE")
writeRaster(rast(threat_stars),paste0(out_directory, threat, "_WIT_raster.tif"))

######  UCI

threat <-  "MagStressInferUCI"

workingsf <- wetland_stress %>% 
  select(UID, RecID, "threat" = all_of(threat))

#testing section
max(workingsf$threat, na.rm=T) 
min(workingsf$threat, na.rm=T) 

#calc section
workingsf = workingsf %>%
  rowwise() %>%
  mutate(threat = abs(min(threat, 0)))

#testing section
max(workingsf$threat, na.rm=T) 
min(workingsf$threat, na.rm=T) 

#calc section
maxthreat <- max(workingsf$threat, na.rm=T)
workingsf = workingsf %>%
  mutate(threat = round(threat/maxthreat,6))   #I am rounding to 6 decimal places here (i think)

#rasterise
in_directory = "D:/Wetland Indices Pilot/Data/"
out_directory <-"D:/Wetland Indices Pilot/Data/"

#read in the wetlands tif as a template
raster_template <- read_stars(paste0(in_directory, "wetlands_WIT_raster.tif"))
raster_template[[1]] <- NA  #need to do this as otherwise end up with an imprint of template raster values

layer2rasterise <- workingsf %>%
  select(threat) 

threat_stars <- st_rasterize(layer2rasterise, template=raster_template, options = "ALL_TOUCHED=TRUE")
writeRaster(rast(threat_stars),paste0(out_directory, threat, "_WIT_raster.tif"))

#end of section that standardises the stress layers

##############################################################################








