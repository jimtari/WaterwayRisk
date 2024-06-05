



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
