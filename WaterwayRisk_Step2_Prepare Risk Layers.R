### This file prepares the species-specific risk layers for integrated risk assessment for a SINGLE asset type
### It requires:
### 1. A threat layer corresponding to a specific threat type
### 2. A species list with species information, including HDM locations, thresholds, & group IDs that can be used to look up relevant Threat-Impact curve shapes
### 3. An associations file giving the T-I curve parameters for relevant taxon groups.
### All input rasters (HDMs, threat and mask layers) are assumed to have consistent extents & projections, corresponding to the 'rasterTemplate' specified

##load packages 

library(tidyverse)
library(sf)
library(stars)
library(terra)
library(mapview)
library(prodlim)

##load custom functions

setwd("O:/SAN_Projects/WaterwaysValuesThreats/WaterwayRisk")

source("Functions for waterway risk assessment calculations.R")

### universal inputs (species list and associations file): set once for all Threats

splist.full=read.csv("O:/SAN_Projects/WaterwaysValuesThreats/WaterwayRisk/WaterwaySpecies_HDMs_Thresholds.csv")
#subset species list to relevant species only ('Include' column by default - but could be based on taxon groups or asset types)
full_list=splist.full%>%subset(Include==1)

#get Associations look up table

associations=read.csv("O:/SAN_Projects/WaterwaysValuesThreats/WaterwayRisk/WaterwaySpecies_ThreatAssociations.csv")


##where should risk layers be saved?
out_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/"
dir.create(file.path(out_directory,"ConditionLayers"))

##where Zonation run files will be saved
Zonation_runDir<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation"

##where Zonation outputs will be saved

Zonation_outDir <-paste(Zonation_runDir,"/Outputs") 

rasterTemplate<-rast("O:/SAN_Projects/WaterwaysValuesThreats/Masks/WaterbodyMask_incEnchee_VicGrid2020.tif")


### Threat-Specific Inputs: change for each Threat and Asset type.

asset.type="Wetland" 

threat_assoc_name<-"Reduced inundation hydrological stress"  #name of Threat in associations file

theat_shortName<-"Reduced_Wetland_Inundation"

#name of threat layer, corresponding to a Threat in the associations file
threat_in <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/HydroStress_wetland_VG20.tif"  

#Set analysis mask layer (could make an asset-specific mask from the threat layer)
analysisMaskFile<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/Waterbodies_VG20.tif"


##### R WORKFLOW STARTS HERE

#make a list variants of the asset.type name (for lookups in other files)
asset.type.alts=c(asset.type,tolower(asset.type))
if(length(grep("s$",asset.type.alts))==0){
  asset.type.alts=c(asset.type.alts,paste0(asset.type.alts,"s"))
}else{
  asset.type.alts=c(asset.type.alts,gsub("s$","",asset.type.alts))  
  }


### get threat layer

threat_raw <- terra::rast(threat_in)
compareGeom(threat_raw,rasterTemplate)

# get waterbody mask
wbody_rast <- terra::rast(analysisMaskFile)
compareGeom(wbody_rast,rasterTemplate)

##Threat layer to Impact layers (via impact curve):
#Create generic 'Impact' maps: one for each curve-shape. 

assocs=associations%>%subset(Threat==threat_assoc_name & Asset%in%asset.type.alts)

curveshapes=unique(assocs[c("M","Y","mu","beta")])

assocs$CurveNumber=row.match(assocs[c("M","Y","mu","beta")],curveshapes)

#Loop thru curve shapes & save Impact and Refuge layers
for(cnum in 1:nrow(curveshapes)){
threat_as_cond<-app(threat_raw,threatToCond,mu.val=curveshapes$mu[cnum],
                                            beta.val=curveshapes$beta[cnum],
                                            M.val=curveshapes$M[cnum],
                                            Y.val=curveshapes$Y[cnum])#%>%classify(rcl=matrix(c(NA,0),nrow=1)) 
##save Impact layer
writeRaster(threat_as_cond,file=file.path(out_directory,"ConditionLayers",paste0(theat_shortName,"_Impact_",cnum,".tif")),overwrite=T)
##save residual value (refugia) condition layer
writeRaster(1-threat_as_cond,file=file.path(out_directory,"ConditionLayers",paste0(theat_shortName,"_Refuge_",cnum,".tif")),overwrite=T)
}


#Set up Zonation runs

setwd(Zonation_runDir)

runname=theat_shortName # name of analysis run (making same as threat name here)

sfname=paste0(runname,"_settings.z5")  # name of the Zonation settings file

condlayers=list.files(file.path(out_directory,"ConditionLayers"),pattern=paste0(theat_shortName,"_Impact"),full.names = T)

condLU=data.frame(index=1:length(condlayers),norm=0,filename=condlayers)

condlink.Name=paste0(runname,"_Impact_CondLinks.txt")

#write condition link file for Zonation: 
write.table(condLU,file=condlink.Name,names=FALSE,row.names=FALSE,quote=FALSE) #Impact layers

#make a refuge value version
condLU$filename=gsub("_Impact_","_Refuge_",condLU$filename)
write.table(condLU,file=gsub("_Impact_","_Refuge_",condlink.Name),names=FALSE,row.names=FALSE,quote=FALSE)

##Make Feature list file

hdmfiles=full_list$WaterwaysHDM
SppNames=full_list$SpeciesName
SppCodes=full_list$TAXON_ID

group=full_list[,paste0(asset.type,"Group")]
output_group <- unique(group)

groupnums=match(group,output_group)

curveNumber=assocs$CurveNumber[match(group,assocs$Group)] #which curve shape (hence _Impact_ layer) to use for each species

conlink=ifelse(is.na(curveNumber),-1,curveNumber)

Feature_list <- data.frame(weight=1, filename=hdmfiles, group=groupnums, condition=conlink,name=SppNames,tr_out=1)

#remove species with no HDM (or VBA buffer) or no Impact curve for this threat
Feature_list=subset(Feature_list,condition>0 & !is.na(filename))

write.table(Feature_list, file=file.path(Zonation_runDir,paste0(runname,"_feature_list.txt")), row.names = F)

## WRITE SETTINGS FILEs
#Risk analysis
settingNames=c("feature list file","analysis area mask layer","hierarchic mask layer","post-processing file","single-feature connectivity link file","condition link file","smoothing factor","stop on plateau","external solution")

settings=data.frame(name=settingNames,value=NA)
rownames(settings)=settingNames

settings["analysis area mask layer","value"]=analysisMaskFile
settings["feature list file","value"]=paste0(runname,"_feature_list.txt")
settings["condition link file","value"] = condlink.Name

usesets=which(!is.na(settings$value))

outcol=paste(settings$name[usesets],settings$value[usesets],sep=" = ")
write.table(outcol,file=sfname,col.names=FALSE,row.names=FALSE,quote=FALSE)

#Settings file for Refugia analysis
settings["condition link file","value"] = gsub("_Impact_","_Refuge_",condlink.Name)
outcol=paste(settings$name[usesets],settings$value[usesets],sep=" = ")
write.table(outcol,file=gsub("_settings","_Refugia_settings",sfname),col.names=FALSE,row.names=FALSE,quote=FALSE)


## WRITE BATCH FILE TO RUN RISK and REFUGE ANALYSES


cat("
@setlocal
@PATH=C:\\Program Files (x86)\\Zonation5;%PATH%

z5 -bwgfc --mode=CAZMAX ",
    paste(sfname,paste0("Output/",paste0(runname,"_risk"))),"\n
z5 -bwgfc --mode=CAZMAX ",
    paste(gsub("_settings","_Refugia_settings",sfname),paste0("Output/",paste0(runname,"_refuge")))
,file=paste0(runname,"_Zrun.cmd"))
    

# RUN ZONATION 
# NB THIS R SESSION WILL NOT BE USEABLE WHILE ZONATION IS RUNNING. DOUBLE CLICK ON THE .cmd file created above to run Z outside of R)

system2(paste0(runname,"_Zrun.cmd"))



# END ANALYSIS. GO TO STEP 3 FILE TO MAKE PRETTY OUTPUTS
# ### STOP HERE


# ### Rendundat Code (I think..)
# ### make a feature list for loading in transformed layers (to get performance & rarity-weighty )
# 
# sfname.load=gsub(".z5","_load.z5",sfname)
# 
# tr_dir=file.path(Zonation_outDir,runname,"transformed_layers")
# rankfile=file.path(Zonation_outDir,runname,"rankmap.tif")
# 
# hdmfiles=list.files(tr_dir,full.names=T,pattern=".tif$")
# SppNames <- gsub(".layer.*$","",gsub(".*Spp","Spp",hdmfiles))
# SppCodes=as.numeric(gsub("Spp","",SppNames))
# 
# groupnums=match(full_list$Group[match(SppCodes,full_list$Taxon_ID)],output_group)
# 
# #Feature_list_load <- data.frame(weight=1, filename=hdmfiles, group=groupnums, name=SppNames)
# 
# write.table(Feature_list_load, file="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/feature_list_load.txt", row.names = F)
# 
# ## WRITE SETTINGS FILEs
# 
# settingNames=c("feature list file","analysis area mask layer","hierarchic mask layer","post-processing file","single-feature connectivity link file","condition link file","smoothing factor","stop on plateau","external solution")
# 
# settings$value=NA
# 
# settings["feature list file","value"]="feature_list_load.txt"
# settings["external solution","value"]=rankfile
# 
# usesets=which(!is.na(settings$value))
# 
# outcol=paste(settings$name[usesets],settings$value[usesets],sep=" = ")
# write.table(settings,file=sfname.load,col.names=FALSE,row.names=FALSE,quote=FALSE)
# 



