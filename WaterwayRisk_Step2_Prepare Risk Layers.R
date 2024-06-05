### This file prepares the species-specific risk layers for integrated risk assessment for a SINGLE asset type
### It requires:
### 1. A threat layer corresponding to a specific threat type
### 2. A species list with species information, including HDM locations, thresholds, & group IDs that can be used to look up relevant Threat-Impact curve shapes
### 3. An associations file giving the T-I curve parameters for relevant taxon groups.
### The names and locations of all inputs are specified in either the Universial settings file, or the Threat-specific settings file, which you must specify below.

##load packages 

library(tidyverse)
library(sf)
library(stars)
library(terra)
library(mapview)
library(prodlim)

## Impact curve function
impactcurve<-function(tval,mu=0.5,beta=0.1,M=1,Y=3.5,B=1/(1+exp(mu/beta)),C=(1+exp(-(1-mu)/beta))/(1-B*(1+exp(-(1-mu)/beta)))){
  ifelse(tval>=Y, M, M*C*(1/(1+exp(-(tval/Y-mu)/beta))-B))
}
## function to apply impact curve function to a raster
threatToCond<-function(inval,mu.val=0.5,beta.val=0.1,M.val=1,Y.val=1){
  outval=impactcurve(tval=inval,mu=mu.val,beta=beta.val,M=M.val,Y=Y.val)
}


#Set working directory
setwd("O:/SAN_Projects/WaterwaysValuesThreats/WaterwayRisk")

##Load universal settings file (check this file)
source("WaterwayRisk_UniversalSettingsFile_TestRuns.R")

##Load Threat-specific settings file (change this file for each Threat and/or variant)
source("WaterwayRisk_ThreatSettings_WetlandHydro_TestRun.R")

##### R WORKFLOW STARTS HERE

#make a list variants of the asset.type name (for lookups in other files)
asset.type.alts=c(asset.type,tolower(asset.type))
if(length(grep("s$",asset.type.alts))==0){
  asset.type.alts=c(asset.type.alts,paste0(asset.type.alts,"s"))
}else{
  asset.type.alts=c(asset.type.alts,gsub("s$","",asset.type.alts))  
  }

## get waterbody mask
wbody_rast <- terra::rast(analysisMaskFile)

## get threat layer
threat_raw <- terra::rast(threat_in)

#check raster match
compareGeom(threat_raw,wbody_rast)

##IF ERROR ABOVE, make sure (in this order)
# 1) the raster specified in 'analysisMaskFile' matches the projection, extent, resolution (and is snapped to):
# "O:/SAN_Projects/WaterwaysValuesThreats/Asset_layers/Streams_Wetlands_Estuaries_categorical_mask_vg20.tif"
#  which should match the HDMs etc
# 2) the raster specific in 'threat_in' matches the raster specified in 'analysisMaskFile'
# threat_raw<- project(threat_raw,wbody_rast)

##Create Impact layers from the Threat layer, via impact curves:
##The threat layer should have positive values only, with zero indicating no threat.
##Make sure the T-I curve parameters correspond to the Threat values in the specified threat layer.

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
write.table(condLU,file=condlink.Name,col.names=FALSE,row.names=FALSE,quote=3) #Impact layers

#make a refuge value version
condLU$filename=gsub("_Impact_","_Refuge_",condLU$filename)
write.table(condLU,file=gsub("_Impact_","_Refuge_",condlink.Name),col.names=FALSE,row.names=FALSE,quote=3)

##Make Feature list file

hdmfiles=file.path(normalHDMs_Dir,basename(full_list$SelectedRaster))

SppNames=full_list$SpeciesName
threshold=ifelse(is.na(full_list$Threshold),0,full_list$Threshold)
SppCodes=full_list$TAXON_ID

group=full_list[,paste0(asset.type,"Group")]
output_group <- unique(group)

groupnums=match(group,output_group)

curveNumber=assocs$CurveNumber[match(group,assocs$Group)] #which curve shape (hence _Impact_ layer) to use for each species

conlink=ifelse(is.na(curveNumber),-1,curveNumber)

Feature_list <- data.frame(weight=1, group=groupnums, condition=conlink,threshold=threshold, name=make.unique(SppNames), filename=hdmfiles,tr_out=saveTransformed)

#remove species with no HDM (or VBA buffer) or no Impact curve for this threat
Feature_list=subset(Feature_list,condition>0 & !is.na(filename))

#check that normalized HDM's exist..

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
## using Zonation flags: a = use analysis area mask, b=pre-processing (rwr, richness),w=use weights (but all 1 here)
## g=use output groups, f=output tansformed layers (=risk layers for each species),
## c= use Condition (i.e. multiply HDMs by Impact layers, corresponding to relevant TI curve shape and Threat layer)
## Adding CAZ2 runs to batch file,but commented out


runname2=paste0(runname,"_CAZ2")

cat("
@setlocal
@PATH=C:\\Program Files (x86)\\Zonation5;%PATH%

z5 -abwgfc --mode=CAZMAX ",
    paste(sfname,paste0("Output/",paste0(runname,"_risk"))),"\n
z5 -abwgfc --mode=CAZMAX ",
    paste(gsub("_settings","_Refugia_settings",sfname),paste0("Output/",paste0(runname,"_refuge"))),"\n
#z5 -abwgfc --mode=CAZ2 ",
    paste(sfname,paste0("Output/",paste0(runname2,"_risk"))),"\n
#z5 -abwgfc --mode=CAZ2 ",
    paste(gsub("_settings","_Refugia_settings",sfname),paste0("Output/",paste0(runname2,"_refuge")))
    
    ,file=paste0(runname,"_Zrun.cmd"))   

##Make a CAZMAX2 version (run serately if needed)
runname2=paste0(runname,"_CAZ2")
cat("
@setlocal
@PATH=C:\\Program Files (x86)\\Zonation5;%PATH%

z5 -abwgfc --mode=CAZ2 ",
    paste(sfname,paste0("Output/",paste0(runname2,"_risk"))),"\n
z5 -abwgfc --mode=CAZ2 ",
    paste(gsub("_settings","_Refugia_settings",sfname),paste0("Output/",paste0(runname2,"_refuge")))
    ,file=paste0(runname,"_Zrun.cmd"))
    
## CHECK THAT ALL HDMs IN FEATURE LIST EXIST: ZONATION WILL NOT RUN IF ANY MISSING
missing.hdms=Feature_list$filename[!file.exists(Feature_list$filename)]

if(length(missing.hdms)>0){
  print("Missing HDMs")
  print(missing.hdms)
  warning("The Normalized HDMs listed above do not exist. Zonation will not run.")
}



# RUN ZONATION 
# NB THIS R SESSION WILL NOT BE USEABLE WHILE ZONATION IS RUNNING. DOUBLE CLICK ON THE .cmd file created above to run Z outside of R)
# Make sure Sei (or wherever running) is not already running Zonation or other memory hungry apps. 


#system2(paste0(runname,"_Zrun.cmd"))   #Runs Zonation in R 



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



