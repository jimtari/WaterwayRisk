### This is a version of the "WaterwatRisk_WaterwayRisk_Step2_Prepare Risk Layers.R" sript, modified to produce the Biodiversity Index
### No threat layer, threat-impact curves or condition layers are used / created in this script

##load packages 

library(tidyverse)
library(sf)
library(stars)
library(terra)
library(mapview)
library(prodlim)




#Set working directory
setwd("O:/SAN_Projects/WaterwaysValuesThreats/WaterwayRisk")

##Load universal settings file (check this file)
source("WaterwayRisk_UniversalSettingsFile_TestRuns.R")

#Set analysis mask layer
analysisMaskFile<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/Waterbodies_VG20.tif"

runname="BiodiversityIndex_Test"  # base name for outputs

##### R WORKFLOW STARTS HERE

# get waterbody mask

#Set up Zonation runs

setwd(Zonation_runDir)


sfname=paste0(runname,"_settings.z5")  # name of the Zonation settings file


##Make Feature list file

hdmfiles=full_list$SelectedRaster
SppNames=full_list$SpeciesName
threshold=ifelse(full_list$HDMorVBA=="HDM",ifelse(is.na(full_list$Threshold),0,full_list$Threshold),0)
SppCodes=full_list$TAXON_ID

group=full_list$TaxonGroup
output_group <- unique(group)

groupnums=match(group,output_group)


Feature_list <- data.frame(weight=1, group=groupnums, threshold=threshold, name=make.unique(SppNames), filename=hdmfiles,tr_out=1)

#remove species with no HDM (or VBA buffer) 
Feature_list=subset(Feature_list,!is.na(filename))

write.table(Feature_list, file=file.path(Zonation_runDir,paste0(runname,"_feature_list.txt")), row.names = F)

## WRITE SETTINGS FILEs
#Risk analysis
settingNames=c("feature list file","analysis area mask layer","hierarchic mask layer","post-processing file","single-feature connectivity link file","condition link file","smoothing factor","stop on plateau","external solution")

settings=data.frame(name=settingNames,value=NA)
rownames(settings)=settingNames

settings["analysis area mask layer","value"]=analysisMaskFile
settings["feature list file","value"]=paste0(runname,"_feature_list.txt")
#settings["condition link file","value"] = condlink.Name

usesets=which(!is.na(settings$value))

outcol=paste(settings$name[usesets],settings$value[usesets],sep=" = ")
write.table(outcol,file=sfname,col.names=FALSE,row.names=FALSE,quote=FALSE)


## WRITE BATCH FILE TO RUN RISK and REFUGE ANALYSES
## using Zonation flags: a = use analysis area mask, b=pre-processing (rwr, richness),w=use weights (but all 1 here)
## g=use output groups, f=output tansformed layers (=risk layers for each species),
## 
## Adding CAZ2 runs to batch file,but commented out


runname2=paste0(runname,"_CAZ2")


cat("
@setlocal
@PATH=C:\\Program Files (x86)\\Zonation5;%PATH%

z5 -abwg --mode=CAZMAX ",
    paste(sfname,paste0("Output/",runname)),"\n
#z5 -abwg --mode=CAZ2 ",
    paste(sfname,paste0("Output/",runname2))
,file=paste0(runname,"_Zrun.cmd")
)
    

# RUN ZONATION 
# NB THIS R SESSION WILL NOT BE USEABLE WHILE ZONATION IS RUNNING. DOUBLE CLICK ON THE .cmd file created above to run Z outside of R)
# Make sure Sei (or wherever running) is not already running Zonation or other memory hungry apps. 

#system2(paste0(runname,"_Zrun.cmd")) #Runs Zonation in R 



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



