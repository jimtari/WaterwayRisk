###############################################################################################
### Threat-specific inputs: set for each Threat (and each Run)                           ######
###############################################################################################

asset.type="Wetland" 

threat_assoc_name<-"Reduced inundation hydrological stress"  #name of Threat in associations file

theat_shortName<-"Reduced_Wetland_Inundation"

#name of threat layer, corresponding to a Threat in the associations file
threat_in <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/Waterbodies_VG20.tif" 
#NB making 'threat' = 1 in all assets: so that result is actually ranking biodiversity values (because threat layer becomes condition layer in Zonation)


#Set analysis mask layer
analysisMaskFile<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/Waterbodies_VG20.tif"

