###############################################################################################
### Universal inputs: set once for all Threats                                           ######
### Loads master species list file, associations/Threat-Impact lookups & directories     ######
### Defines T-I curve Functions                                                          ######     
###############################################################################################

library(readxl)

splist.full=read_excel("O:/SAN_Projects/WaterwaysValuesThreats/WaterwaySpecies_HDMs_Thresholds.xlsx",sheet=1,na="NA")

#subset species list to relevant species only ('Include' column by default - but could be based on taxon groups or asset types)
full_list=splist.full%>%subset(Include==1)

#get Associations look up table

associations=read_excel("O:/SAN_Projects/WaterwaysValuesThreats/WaterwaySpecies_HDMs_Thresholds.xlsx",sheet=2,na="NA")

##where should risk layers be saved?
out_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/"
dir.create(file.path(out_directory,"ConditionLayers"))

##where Zonation run files will be saved
Zonation_runDir<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation"

##where Zonation outputs will be saved

Zonation_outDir <-paste(Zonation_runDir,"/Outputs") 





