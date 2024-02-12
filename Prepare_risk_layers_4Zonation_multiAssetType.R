### This file prepares the risk (threat*biodiversity) layers for Zonation.
## and writes the necessary files. This code uses threat layers prepared by 
## "Prepare_threats_4zonation_WIT.R"



library(tidyverse)
library(sf)
library(stars)
library(terra)
library(mapview)
library(doParallel)

#impact curve function


impactcurve<-function(tval,mu=0.5,beta=0.1,M=1,Y=3.5,B=1/(1+exp(mu/beta)),C=(1+exp(-(1-mu)/beta))/(1-B*(1+exp(-(1-mu)/beta)))){
  ifelse(tval>=Y, M, M*C*(1/(1+exp(-(tval/Y-mu)/beta))-B))
}

#function to apply impact curve function to a raster

threatToCond<-function(inval,type=2,mu.val=0.5,beta.val=0.1,M.val=1,Y.val=3.5){
  if(type==1){mu.val=0;beta.val=0.1;M.val=1;Y.val=3.5}
  if(type==3){mu.val=1;beta.val=1;M.val=1;Y.val=3.5}
  if(type==4){mu.val=1;beta.val=0.1;M.val=1;Y.val=3.5}
  impactcurve(tval=inval,mu=mu.val,beta=beta.val,M=M.val,Y=Y.val)
}
  


#absolute of negative stress values only

absStress<-function(stress,posneg="neg"){   #posneg='neg' to make positive stress zero, posneg='pos' to make negative stress zero
  stressout=abs(stress)
  if(posneg=="neg"){stressout=stressout*(stress<0)}
  if(posneg=="pos"){stressout=stressout*(stress>0)}
  stressout
}

#########################################################################################
## Define output groups
########################################################################################


 output_group <- c("Plants", "Birds", "Amphibians","Reptiles","Fish","Inverts","Other")  #Options are: "Herps", "Plants", "Birds", "Aquatic", c("Herps", "Plants", "Birds", "Aquatic")
 output_label <- "All"

# output_group <- "Herps"
# output_label <- "Herps"
 

 


#########################################################################################
## Read in the threat layer(s) 
########################################################################################

threat_in <-"HydroStress_dry"
threat_out<- threat_in
threat_in_wetland <- paste0(threat_in,"_raster_wetland.tif")
threat_in_reach <- paste0(threat_in,"_raster_reach.tif")
threat_in_estuary <- paste0(threat_in,"_raster_estuary.tif")  


in_directory<-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/"
out_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/"
dir.create(file.path(out_directory,"ConditionLayers"))
dir.create(file.path(out_directory,"HDMs"))

wetland_threat_raw <- terra::rast(file.path(in_directory, threat_in_wetland))
reach_threat_raw <- terra::rast(file.path(in_directory, threat_in_reach))
estuary_threat_raw <- terra::rast(file.path(in_directory, threat_in_estuary))

wbody_rast <- terra::rast(paste0(in_directory, "Waterbodies_raster.tif"))

#stress to condition (via impact curve):
wetland_threat_scaled<-app(wetland_threat_raw,absStress)
reach_threat_scaled<-4*(10-reach_threat_raw)/10 

writeRaster(wetland_threat_scaled,file=file.path(out_directory,paste0(threat_out,"_wetland.tif")),overwrite=T)
writeRaster(reach_threat_scaled,file=file.path(out_directory,paste0(threat_out,"_reach.tif")),overwrite=T)

##combinations of impacts (asset types X curve types)

for(wettype in c(2,4)){
  for(reachtype in c(2,4)){
    esttype=2

#weltands
#ignore positive stress, get absolute value of negative stress
wetland_cond<-app(wetland_threat_scaled,threatToCond,type=wettype)%>%classify(rcl=matrix(c(NA,0),nrow=1))

#reaches
 #scaling FSR stress to 0 (no stress) to 7 (max stress), from 10 no stress to 0 max stress
reach_cond<-app(reach_threat_scaled,threatToCond,type=reachtype)%>%classify(rcl=matrix(c(NA,0),nrow=1))

#estuary
estuary_threat_scaled=estuary_threat_raw%>%classify(rcl=matrix(c(NA,0),nrow=1)) # no drying stress for estuary
estuary_cond=estuary_threat_scaled

threat_as_cond=app(c(wetland_cond,reach_cond,estuary_cond),max,na.rm=T)  #combined condition (threat impact) is max of asset type conditions

writeRaster(threat_as_cond,file=file.path(out_directory,"ConditionLayers",paste0(threat_out,"_Cond_",wettype,reachtype,esttype,".tif")),overwrite=T)

}}

#########################################################################################
## Clip all the HDMs selected earlier (in filenames) to the wetlands, 
########################################################################################



library(terra)
library(doParallel)

wbody_rast=rast("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/Waterbodies_raster.tif")>0

feature_in_directory <- "O:/SAN_Projects/Wetland_HDMs/cropped_and_maskedOriginal_models/MaskedWetlandSppBig"
filenames=list.files(feature_in_directory,pattern=".tif$")
outnames=gsub("Masked_","",filenames)
hdm_out_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/HDMs"

paths=.libPaths()

cl <- makeCluster(20,outfile="")
registerDoParallel(cl)
hdmsums<-foreach(spp=iter(filenames[duds]),.packages=c("terra"),.combine="rbind",.errorhandling="remove")%dopar%{
  .libPaths(paths)
  wbody_rast=rast("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/Waterbodies_raster.tif")>0
  filename <- file.path(feature_in_directory,spp)
  outname<-gsub("Masked_","",spp)
  feat_raw <- terra::rast(filename)
  rawtot <- global(feat_raw, 'sum', na.rm=T)
  feat <- project(feat_raw,wbody_rast, method="bilinear")*(wbody_rast>0)
  tot <- global(feat, 'sum', na.rm=T) #find the sum value of all cells in the raster

  terra::writeRaster(round(feat),file.path(hdm_out_directory,outname), NAflag=0,gdal="COMPRESS=DEFLATE",overwrite=TRUE,datatype="INT1U")
  data.frame(taxonID=gsub(".tif$","",gsub(".*Spp","",filename)),hdm=outname,HDMsum_wwy=tot[1,1],HDMsum=rawtot[1,1])

}

stopCluster(cl)


duds=which(!file.exists(file.path(hdm_out_directory,gsub("Masked_","",filenames))))
if(length(duds)>0){hdmsums1=hdmsums}  # now redo above code using iter(filenames[duds]), and maybe errorhandling='stop' to debug

## after rerunning duds (& when lenth(duds)==0)

hdmsums=rbind(hdmsums1,hdmsums)

write.csv(hdmsums,file=file.path(hdm_out_directory,"HDMsums.csv"),row.names=F)

##test hdms
hdmfiles=list.files(hdm_out_directory,full.names=T)
rr=rast(hdmfiles)

dudhdms=which(is.na(apply(minmax(rr),2,diff)))

hdmfiles=hdmfiles[-dudhdms]



#########################################################################################
## Make Z5 feature files that includes groups info as well as location and name of files 
##########################################################################################


threat_in <-"HydroStress_dry"

sfname="settings_vic_risk_test.z5"

zfiles_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation"

hdm_out_directory <-"O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/HDMs"

hdmfiles=list.files(hdm_out_directory,full.names=T)
shortnames <- gsub(".tif$","",gsub(".*Spp","Spp",hdmfiles))
codes=as.numeric(gsub("Spp","",shortnames))
SppNames=gsub(".tif$","",gsub(".*Spp","Spp",hdmfiles))


condlayers=list.files("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/ConditionLayers",pattern=threat_in,full.names = T)

condLU=data.frame(index=1:length(condlayers),norm=0,filename=condlayers)

##Make Feature list file

hdmfiles=list.files(hdm_out_directory,full.names=T)
SppNames=gsub(".tif$","",gsub(".*Spp","Spp",hdmfiles))
SppCodes=as.numeric(gsub("Spp","",SppNames))

full_list <- read.csv("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/HDM_models_species_Appendix1plus.csv")

full_list$Group=gsub("_.*","",full_list$AVIRA_group)
output_group <- unique(full_list$Group)

groupnums=match(full_list$Group[match(SppCodes,full_list$Taxon_ID)],output_group)

associations=full_list$hydro_stress_assoc[match(SppCodes,full_list$Taxon_ID)]
conlink=ifelse(is.na(associations),-1,ifelse(associations=="med",4,1))

Feature_list <- data.frame(weight=1, filename=hdmfiles, group=groupnums, condition=conlink,name=SppNames,tr_out=1)

write.table(Feature_list, file="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/feature_list.txt", row.names = F)

## WRITE SETTINGS FILEs
settingNames=c("feature list file","analysis area mask layer","hierarchic mask layer","post-processing file","single-feature connectivity link file","condition link file","smoothing factor","stop on plateau","external solution")

settings=data.frame(name=settingNames,value=NA)
rownames(settings)=settingNames

settings["feature list file","value"]="feature_list.txt"
settings["condition link file"] = "HydroStress_dry_CondLinks.txt"

usesets=which(!is.na(settings$value))

outcol=paste(settings$name[usesets],settings$value[usesets],sep=" = ")
write.table(settings,file=sfname,col.names=FALSE,row.names=FALSE,quote=FALSE)




#### make a feature list for loading in transformed layers (to get performance)


runname="vic_test_q10"
sfname.load=gsub(".z5","_load.z5",sfname)

tr_dir=file.path(zfiles_directory,"Outputs",runname,"transformed_layers")
rankfile=file.path(zfiles_directory,"Outputs",runname,"rankmap.tif")

hdmfiles=list.files(tr_dir,full.names=T,pattern=".tif$")
SppNames <- gsub(".layer.*$","",gsub(".*Spp","Spp",hdmfiles))
SppCodes=as.numeric(gsub("Spp","",SppNames))

groupnums=match(full_list$Group[match(SppCodes,full_list$Taxon_ID)],output_group)

Feature_list_load <- data.frame(weight=1, filename=hdmfiles, group=groupnums, name=SppNames)

write.table(Feature_list_load, file="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/feature_list_load.txt", row.names = F)

## WRITE SETTINGS FILEs

settingNames=c("feature list file","analysis area mask layer","hierarchic mask layer","post-processing file","single-feature connectivity link file","condition link file","smoothing factor","stop on plateau","external solution")

settings$value=NA

settings["feature list file","value"]="feature_list_load.txt"
settings["external solution","value"]=rankfile

usesets=which(!is.na(settings$value))

outcol=paste(settings$name[usesets],settings$value[usesets],sep=" = ")
write.table(settings,file=sfname.load,col.names=FALSE,row.names=FALSE,quote=FALSE)





#########################################################
######   SUMMARIZE RESULTS  #############################
##This is all experimental / mess. Tidy up.
#################
#extract results / inputs to shapefile...

isc=st_read("O:/SAN_Projects/WaterwaysValuesThreats/Asset_layers/ISC2010_FSR2020_stats.shp")
isc2=st_simplify(isc,dTolerance = 100)

wetland <- readRDS(paste0(getwd(), "/PilotV4/wetland_stress.rds"))

iscbuff<-st_buffer(isc,dist=75)

iscbuff$threat<-(iscbuff$Q90s+iscbuff$PZDs)/2

iscbuff$threat<-iscbuff$Q10s# +iscbuff$PZDs)/2

iscbuff$impact=threatToCond(4*(10-iscbuff$threat)/10)

egsp=rast("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Data/HDMs/River_Blackfish_Spp4949.tif")

sppn=egsp/sum(values(egsp),na.rm=T)

sppv=terra::extract(sppn,st_transform(iscbuff,crs(sppn)),fun="sum",na.rm=T)

iscbuff$Spp4949=sppv$Masked_River_Blackfish_Spp4949

iscbuff$River_blackfish_Spp4949=iscbuff$Spp4949/sum(iscbuff$Spp4949,na.rm=T)

iscbuff$Spp4949_risk=iscbuff$Spp4949*iscbuff$impact

iscbuff$Risk_Spp4949=iscbuff$River_blackfish_Spp4949*iscbuff$impact



#curves=read.table("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/Outputs/vic_test2_isc/feature_curves.csv",header=T)

#isz=which(curves[nrow(curves),-1]==0)
#rstack=rast(zlist$filename)
#rbin=rstack>0
#nspp=app(rstack>0,sum,na.rm=T)

rwr=rast("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/Outputs/vic_test_q10/pre_processing/rsr.pos.tif")
zrank=rast("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/Outputs/vic_test_q10/rankmap.tif")

rwrvals=values(rwr)
zvals=values(zrank)

haveboth=which(!is.na(rwrvals) & !is.na(zvals) & rwrvals>0)

rwrrank=rank(rwrvals,na.last = 'keep',ties.method = 'random')
rwrrank=rwrrank/max(rwrrank,na.rm=T)
rs=sample(haveboth,100000)

gdat=data.frame(rwr=rwrvals[rs],zrank=zvals[rs],rwrrank=rwrrank[rs])




library(ggplot2)

p1=ggplot(gdat)+geom_point(aes(x=zrank,y=rwr))#+scale_y_continuous(trans='log2')
p2=ggplot(gdat)+geom_point(aes(x=zrank,y=rwrrank))

p1+theme_bw()+ylab("Sum of Species Risks")+xlab("Risk Rank (Zonation)")


plot(rwrvals[rs]~zvals[rs])

plot(rwrrank[rs]~zvals[rs])

Risk_total=terra::extract(rwr,st_transform(iscbuff,crs(rwr)),fun="sum",na.rm=T,touches=T)$rsr.pos
Risk_mean=terra::extract(rwr,st_transform(iscbuff,crs(rwr)),fun="mean",na.rm=T,touches=T)$rsr.pos
Risk_Zrank=terra::extract(zrank,st_transform(iscbuff,crs(zrank)),fun="mean",na.rm=T,touches=T)$rankmap
Risk_Zrank_max=terra::extract(zrank,st_transform(iscbuff,crs(zrank)),fun="max",na.rm=T,touches=T)$rankmap
Risk_Zrank_sum=terra::extract(zrank,st_transform(iscbuff,crs(zrank)),fun="sum",na.rm=T,touches=T)$rankmap

wetland.t=st_transform(wetland,crs(rwr))

Risk_total=terra::extract(rwr,wetland.t,fun="sum",na.rm=T,touches=T)$rsr.pos
Risk_mean=terra::extract(rwr,wetland.t,fun="mean",na.rm=T,touches=T)$rsr.pos
Risk_Zrank=terra::extract(zrank,wetland.t,fun="mean",na.rm=T,touches=T)$rankmap
Risk_Zrank_max=terra::extract(zrank,wetland.t,fun="max",na.rm=T,touches=T)$rankmap
Risk_Zrank_sum=terra::extract(zrank,wetland.t,fun="sum",na.rm=T,touches=T)$rankmap


wetl_simp=st_simplify(wetland,dTolerance=100)
wetl_cent=st_centroid(wetland)

wetl_cent$Risk_All=Risk_total
wetl_cent$Risk_rank=Risk_Zrank
wetl_cent$Risk_mean=Risk_mean
wetl_cent$Risk_Zmax=Risk_Zrank_max
wetl_cent$Risk_Zsum=Risk_Zrank_sum
wetl_cent$MagStr=wetl_cent$MagStressInfer
wetl_cent=wetl_cent[nchar(names(wetl_cent))<11]

wetl_gps <- wetl_cent %>% group_by(IWG_ID) %>% summarise() %>% st_centroid()

wgtmean=function(x,w){sum(x*w)/sum(w)}

gpstats<-st_drop_geometry(wetl_cent)%>%group_by(IWG_ID)%>%summarise(across(starts_with("Risk"),.fns=~weighted.mean(.x,HECTARES,na.rm=T)))%>%data.frame()
gpstats$Risk_All=tapply(wetl_cent$Risk_All,wetl_cent$IWG_ID,sum,na.rm=T)

wetl_gps=cbind(wetl_gps,gpstats)

wetl_gps$AREA=tapply(wetl_cent$HECTARES,wetl_cent$IWG_ID,sum,na.rm=T)

st_write(wetl_cent,dsn="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/RiskDemo_wetland_points.shp",append=FALSE)
st_write(wetl_gps,dsn="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/RiskDemo_wetlandGroup_points.shp",append=FALSE)

riskdemo=iscbuff%>%st_simplify(dTolerance=10)

riskdemo2=st_buffer(riskdemo,dist=1000)

riskdemo2$threat=(10-riskdemo$threat)

riskdemo2$Risk_All=Risk_total
riskdemo2$Risk_rank=Risk_Zrank
riskdemo2$Risk_mean=Risk_mean
riskdemo2$Risk_Zmax=Risk_Zrank_max
riskdemo2$Risk_Zsum=Risk_Zrank_sum


isc2$Risk_All=Risk_total
isc2$Risk_rank=Risk_Zrank
isc2$Risk_mean=Risk_mean
isc2$Risk_Zmax=Risk_Zrank_max
isc2$Risk_Zsum=Risk_Zrank_sum


riskdemo2=riskdemo2[names(riskdemo2)!="Spp4949_risk"]

names(riskdemo2)=gsub("River_blackfish_Spp4949","HDM_spp",names(riskdemo2))
names(riskdemo2)=gsub("Risk_Spp4949","Risk_spp",names(riskdemo2))


st_write(riskdemo2,dsn="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/RiskDemo_reach.shp",append=FALSE)

st_write(isc2,dsn="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/RiskDemo_reach_lines_q10.shp",append=FALSE)

pal.spect<-function(n){hcl.colors(n,palette='Spectral',rev=T)}
pal.virid<-function(n){hcl.colors(n)}


p.threat=ggplot(riskdemo2)+geom_sf(aes(fill=threat),color=NA)+
  scale_fill_viridis_c(option="B")+coord_sf(datum = NA)+theme_bw()+
  theme(panel.background=element_blank(),panel.border = element_blank(),legend.position=c(0.05,0.75),legend.title=element_blank(),legend.background=element_blank())

p.impact=ggplot(riskdemo2)+geom_sf(aes(fill=impact),color=NA)+
  scale_fill_viridis_c(option="A")+coord_sf(datum = NA)+theme_bw()+
  theme(panel.background=element_blank(),panel.border = element_blank(),legend.position=c(0.05,0.75),legend.title=element_blank(),legend.background=element_blank())

p.hdm=ggplot(riskdemo2)+geom_sf(aes(fill=HDM_spp),color=NA)+
  scale_fill_viridis_c()+coord_sf(datum = NA)+theme_bw()+
  theme(panel.background=element_blank(),panel.border = element_blank(),legend.position=c(0.05,0.75),legend.title=element_blank(),legend.background=element_blank())


p.risk=ggplot(riskdemo2)+geom_sf(aes(fill=Risk_spp),color=NA)+
  scale_fill_viridis_c(option="H")+coord_sf(datum = NA)+theme_bw()+
  theme(panel.background=element_blank(),panel.border = element_blank(),legend.position=c(0.05,0.75),legend.title=element_blank(),legend.background=element_blank())

p.risk_all=ggplot(riskdemo2)+geom_sf(aes(fill=Risk_All),color=NA)+
  scale_fill_viridis_c(option="H")+coord_sf(datum = NA)+theme_bw()+
  theme(panel.background=element_blank(),panel.border = element_blank(),legend.position=c(0.05,0.75),legend.title=element_blank(),legend.background=element_blank())


p.risk_rank=ggplot(riskdemo2)+geom_sf(aes(fill=Risk_rank),color=NA)+
  scale_fill_viridis_c(option="H")+coord_sf(datum = NA)+theme_bw()+
  theme(panel.background=element_blank(),panel.border = element_blank(),legend.position=c(0.05,0.75),legend.title=element_blank(),legend.background=element_blank())

p.risk_mean=ggplot(riskdemo2)+geom_sf(aes(fill=Risk_mean),color=NA)+
  scale_fill_viridis_c(option="H")+coord_sf(datum = NA)+theme_bw()+
  theme(panel.background=element_blank(),panel.border = element_blank(),legend.position=c(0.05,0.75),legend.title=element_blank(),legend.background=element_blank())



library(ggpubr)
ggexport(p.risk_all,res=250,width=1000,height=1000,filename="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Risk_All.png",pointsize=5)
ggexport(p.risk_rank,res=250,width=1000,height=1000,filename="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Risk_Rank.png",pointsize=5)
ggexport(p.risk_mean,res=250,width=1000,height=1000,filename="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Risk_Mean.png",pointsize=5)
ggexport(p.risk,res=250,width=1000,height=1000,filename="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Risk.png",pointsize=5)
ggexport(p.threat,res=250,width=1000,height=1000,filename="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Threat.png",pointsize=5)
ggexport(p.impact,res=250,width=1000,height=1000,filename="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Impact.png",pointsize=5)
ggexport(p.hdm,res=250,width=1000,height=1000,filename="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/HDM.png",pointsize=5)




####
iscbuff$Risk_All=isc2$Risk_All
sumR.isc=rasterize(st_transform(iscbuff,crs(rwr)),rwr,field="Risk_All")
wetland.t$Risk_All=wetl_cent$Risk_All
sumR.wet=rasterize(wetland.t,rwr,field="Risk_All")
sumR=app(c(sumR.isc,sumR.wet),max,na.rm=T)

writeRaster(sumR,file="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/Outputs/vic_test_q10/sumRisk_wbody.tif")

iscbuff$Risk_rank=isc2$Risk_rank
Zmean.isc=rasterize(st_transform(iscbuff,crs(rwr)),rwr,field="Risk_rank")
wetland.t$Risk_rank=wetl_cent$Risk_rank
Zmean.wet=rasterize(wetland.t,rwr,field="Risk_rank")
Zmean=app(c(Zmean.isc,Zmean.wet),max,na.rm=T)

writeRaster(Zmean,file="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/Outputs/vic_test_q10/Zrank_wbody.tif")

Rsums=c(iscbuff$Risk_All,wetland.t$Risk_All)
Zranks=c(iscbuff$Risk_rank,wetland.t$Risk_rank)
areas=c(st_area(iscbuff),st_area(wetland.t))

Rsum.area=cumsum(areas[order(Rsums,decreasing=T)])
Zrank.area=cumsum(areas[order(Zranks,decreasing=T)])

Rsum.propA=Rsum.area/max(Rsum.area)
Zrank.propA=Zrank.area/max(Zrank.area)


prop.wbody=c(1:length(areas))/length(areas)

write.csv(data.frame(Rsum=Rsum.propA,Zrank=Zrank.propA,Wbody=prop.wbody),file="O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/prop.waterbodies.csv")


#######################

zlist=read.table("O:/SAN_Projects/WaterwaysValuesThreats/Waterway Risk Pilot/Zonation/feature_list.txt", header=T)
library(doParallel)

paths <- .libPaths()

cfun=function(v1,v2){v1+v2}

cl <- makeCluster(5)
registerDoParallel(cl)

foreach(spp=iter(zlist$filename[1:5]),.packages=c("terra"),.errorhandling="remove")%dopar%{
  .libPaths(paths)
  print(spp)
  #ifelse(is.na(values(r)),0,ifelse(values(r)>0,1,0))
}
stopCluster(cl)  
