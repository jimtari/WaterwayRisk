###############################################################################################
### Universal inputs: set once for all Threats                                           ######
### Loads master species list file, associations/Threat-Impact lookup                    ######
### Sets Directories, and Defines T-I curve Functions                                    ######     
###############################################################################################

library(tidyverse)
library(readxl)


splist.full=read_excel("O:/SAN_Projects/WaterwaysValuesThreats/WaterwaySpecies_HDMs_Thresholds.xlsx",sheet=1,na="NA")%>%data.frame()

#subset species list to relevant species only ('Include' column by default - but could be based on taxon groups or asset types)
full_list=splist.full%>%subset(Include==1)

#get Associations look up table

associations=read_excel("O:/SAN_Projects/WaterwaysValuesThreats/WaterwaySpecies_HDMs_Thresholds.xlsx",sheet=2,na="NA")%>%data.frame()

##where should risk layers be saved?
out_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/"
dir.create(file.path(out_directory,"ConditionLayers"))

##where Zonation run files will be saved
Zonation_runDir<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation"

##where Zonation outputs will be saved

Zonation_outDir <-paste(Zonation_runDir,"/Outputs") 

normalHDMs_Dir<-"E:/Zonation/WaterwaysValuesThreats/HDMs_normalised"

#impact curve function

impactcurve<-function(tval,mu=0.5,beta=0.1,M=1,Y=3.5,B=1/(1+exp(mu/beta)),C=(1+exp(-(1-mu)/beta))/(1-B*(1+exp(-(1-mu)/beta)))){
  ifelse(tval>=Y, M, M*C*(1/(1+exp(-(tval/Y-mu)/beta))-B))
}

#function to apply impact curve function to a raster

threatToCond<-function(inval,mu.val=0.5,beta.val=0.1,M.val=1,Y.val=1){
  outval=impactcurve(tval=inval,mu=mu.val,beta=beta.val,M=M.val,Y=Y.val)
}
