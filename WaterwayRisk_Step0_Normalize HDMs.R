
library(terra)
library(doParallel)

#Set working directory
setwd("O:/SAN_Projects/WaterwaysValuesThreats/WaterwayRisk")

##Load universal settings file (check this file)
source("WaterwayRisk_UniversalSettingsFile_TestRuns.R")

assetMasterFile<-"O:/SAN_Projects/WaterwaysValuesThreats/Asset_layers/Streams_Wetlands_Estuaries_categorical_mask_vg20.tif"

analysisMaskFile<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/Waterbodies_VG20.tif"

asset=rast(assetMasterFile)
maskfile=rast(analysisMaskFile)
compareGeom(asset,maskfile)

type.lu=levels(asset)[[1]]



spcodes=full_list$TAXON_ID[!is.na(full_list$SelectedRaster) & file.exists(full_list$SelectedRaster)]

cl <- makeCluster(20)
registerDoParallel(cl)
HDMsums<-foreach(spp=iter(spcodes),.packages=c("terra","tcltk","dplyr"),.errorhandling="remove",.combine=rbind)%dopar%{

  asset=rast(assetMasterFile)
  maskfile=rast(analysisMaskFile)
  type.lu=levels(asset)[[1]]
  
  #threshold & clip to analysisMaskFile

#rasname=full_list$WaterwaysHDM[match(spp,full_list$TAXON_ID)]
rasname=full_list$SelectedRaster[match(spp,full_list$TAXON_ID)]  
ras=rast(rasname)
thresh=NA
thresh=full_list$Threshold[match(spp,full_list$TAXON_ID)]
if(is.na(thresh) | full_list$HDMorVBA[match(spp,full_list$TAXON_ID)]!="HDM"){thresh=0}

sumRaw=sum(values(ras),na.rm=T)

#threshold
ras=ras*(ras>=thresh)

sumT=sum(values(ras),na.rm=T)

#clip to waterways mask
ras=ras*maskfile

sumClip=sum(values(ras),na.rm=T)


#get HDM sums within each asset class

type.lu$sumHDM=NA

for(rid in 1:nrow(type.lu)){  

  typeid=type.lu$value[rid]
  type.lu$sumHDM[rid]=sum(values(ras*(asset==typeid)),na.rm=T)
}
  tsums=type.lu$sumHDM
  names(tsums)=type.lu$Waterway_Type
  
 

  #normalize
  outras=ras/sumClip
  
  #write to NormalisedHDMfile
  outname=gsub(dirname(rasname),normalHDMs_Dir,rasname)
  writeRaster(outras,file=outname,overwrite=T,gdal="COMPRESS=DEFLATE")
  
  if(!exists("pb0")){pb0 <- tkProgressBar("Clipping Waterway Rasters", min=0, max=length(spcodes))}
  setTkProgressBar(pb0, value=match(spp,spcodes))
  
  #compile sums
  c(TAXON_ID=spp,HDM=rasname,normalizedHDM=outname,sumRaw=sumRaw,sumT=sumT,sumClip=sumClip,tsums)

}  
stopCluster()


write.csv(HDMsums,file="O:/SAN_Projects/WaterwaysValuesThreats/WaterwayRisk/HDMsums_byAssetType.csv")
  