install.packages('writexl')
library(tidyverse)
library(magrittr)
library(stringr)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggforce)
library(corrplot)
library(vegan)

library(readxl)

#########################MALASPINA#######################################
MP_VP <- read.delim("C:/Users/Desktop/TFM/malaspina/MP-RGC-v4-vertProf-length-norm-counts-LuxA.tbl")
MP_rgc_VP <- read.delim("C:/Users/Desktop/TFM/malaspina/MP-RGC-v4-vertProf-length-norm-scg-norm-counts-luxA.tbl")
malaspina.codis <- read_excel("C:/Users/Desktop/TFM/malaspina/media-1.xlsx", sheet = 3)

MP_rgc_VP_t <- t(MP_rgc_VP) #transposar taula pq quedi com malaspina_codis
MP_rgc_VP_t<-as.data.frame(MP_rgc_VP_t)
colnames(MP_rgc_VP_t) <- MP_rgc_VP_t[1,] #posar primera fila com a nom de columna
MP_rgc_VP_t<-MP_rgc_VP_t[-1,] #eliminar primera fila pq ja esta com a nom de columnes
MP_rgc_VP_t<-data.frame(MP_rgc_VP_t)
MP_rgc_VP_t<-rownames_to_column(MP_rgc_VP_t, "MPCode") #passar noms de fila a columna nova que es digui MPCode
MP_abund<-left_join(MP_rgc_VP_t, malaspina.codis, by="MPCode") #fusionar malaspina.codis amb les abundancies de cada gen
MP_abund$depth #per veure el contingut de la columna depth per exemple
MP_abund$sampling_station
rownames(MP_rgc_VP_t) <- MP_rgc_VP_t[,1]
#MP_rgc_VP_t<-MP_rgc_VP_t[,-1]
malaspina.codis<-as.data.frame(malaspina.codis)
rownames(malaspina.codis) <- malaspina.codis[,1]
colnames(malaspina.codis)
rownames(MP_rgc_VP_t)
view(MP_rgc_VP_t)


x<-tidyr::gather(MP_abund, variables, values, MP.RGC.v4_001193594:MP.RGC.v4_067950009)
x$values
x$variables
head(x)
names(x)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
superpalette<-c("#8D4C6A", #cyano-cyano
                "#c2abc9", #gamma-unk
                "#b27300", #bacte-unk
                "#FFDF26", #unk-unk (0.2 no)
                "#626591", #alpha-cladeI
                "#0000ff", #alpha-cladeII (0.2)
                "#d8e6f3", #alpha-cladeIII
                "#66b2b2", #alpha-rhodo
                "#00004c", #alpha-stappia
                "#005b96", #alpha-AEGEAN (0.2)
                "#006666", #alpha-hypho (0.2 no)
                "#377EB8", #alpha-parv (0.2 no)
                "#00ff99", #alpha-sphingo
                "#4DAF4A", #gamma-comamo
                "#74ee15", #gamma-methyl
                "#028900", #gamma-burkhol
                "#740001", #gamma-oxalo  (0.2 no)
                "#ff2032", #gamma-morax
                "#4F507F", #gamma-pseudomona
                "#99cc99", #gamma-litorico
                "#6073b6", #gamma-spongiibact  (0.2 no)
                "#fe9360", #bacter-flavobacteria
                "#ff00b4", #bacter-NS9
                "#f3fc64", #Bacter-Crocinito
                "#b12ef7", #bacter-cryomo
                "#a7dbd8", #bacter-saprospir  (0.2 no)
                "#ea80fc", #bacter-chitinoph
                "#310e59", #bacter-NS11-12
                "#B15A7A", #bacter-cyclobact
                "#ff2e65", #bacter-spirosom (0.2)
                "#e6e8ea", #bacter-hymeno (0.2 no)
                "#FF7F00", #rhodothermia-balneolacae
                "#dab600", #actino-sporichthya
                "#ff80ed", #actino-microbact
                "#F15E75", #verruco-rubrital (0.2 no)
                "#FFFF33", #verruco-punicei
                "#ff0000", #amatimona-amatimona
                "#4169E1", #bacillii-planococc (0.2 no)
                "#0000CD") #z_other
palette<-colorRampPalette(cbbPalette)
palette_2<-colorRampPalette(superpalette)


##Abundàncies totals
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>%
  data.frame()%>%
  #subset(sampling_station=="MH_019") %>%
  na.exclude()%>% #excloure missing values
  select(c("sampling_station", "depth", "layer", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(values), fill = variables)) +
  #geom_point() + geom_line()+
  geom_bar(position ="stack", stat="identity")+#graph amb abund totals (position ="stack")
  scale_fill_manual(values=palette_2(90))+    
  # scale_y_log10() +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(x="Layers", y ="LuxA gene total abundance")

##Abundàncies relatives
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>% 
  data.frame()%>%
  # subset(sampling_station=="MH_019") %>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(values), fill = variables)) +
  #geom_point() + geom_line()+
  geom_bar(position ="fill", stat="identity")+ #graph amb abund relatives (position="fill") 
  scale_fill_manual(values=palette_2(90))+
  # scale_y_log10() +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(x="Layers", y ="LuxA gene relative abundance")

##Gràfic depth/long 
names(MP_abund)
rownames(MP_abund)
x<-tidyr::gather(MP_abund, variables, values, MP.RGC.v4_001193594:MP.RGC.v4_067950009)
x$values
x$variables
head(x)
names(x)

x$layer<-factor(x$layer, levels=c("epipelagic", "DCM", "mesopelagic", "bathypelagic" ))

##Distribució vertical
vertical_palette <- c("epipelagic" = "#FFFF33","#b7db70","#cce5ff","#cccccc")
x %>% 
  data.frame()%>%
  #subset(sampling_station=="MH_019") %>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "longitude", "variables", "values")) %>%
  ggplot(aes(x=longitude, y =depth, size = values, color = sampling_station, shape = layer)) +
  geom_rect(inherit.aes=F, aes(xmin=-Inf, xmax=Inf, ymin=0, max = 200), fill = "#ffff99", alpha =0.02)+
  geom_rect(inherit.aes=F, aes(xmin=-Inf, xmax=Inf, ymin=200, max = 1000), fill = "#cce5ff", alpha =0.02)+
  geom_rect(inherit.aes=F, aes(xmin=-Inf, xmax=Inf, ymin=1000, max = 4000), fill = "#cccccc", alpha =0.02)+
  geom_point(alpha=0.7) +
  scale_shape_manual(values = c("epipelagic"=16,"DCM"= 13, "mesopelagic"= 16, "bathypelagic"=16))+
  #scale_y_reverse()+
  #scale_y_log10()+
  scale_y_continuous(trans=trans_reverser('log10'))+
  scale_color_manual(values=palette_2(11))+
  theme_bw()+
  #theme(aspect.ratio = 4/20)+
  theme(legend.position="none")+
  labs(x="Longitude", y ="Depth (m)")

##Temperatura
#graph parameters
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>% 
  data.frame()%>%
  # subset(sampling_station=="MH_019") %>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "temperature", "variables", "values")) %>%
  ggplot(aes(x=layer, y =temperature, size = values)) +
  geom_point(color = "#FF6633", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Temperature (ºC)", x ="Layer")


##Salinitat
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "salinity", "variables", "values")) %>%
  ggplot(aes(x=layer, y =salinity, size = values)) +
  geom_point(color = "#FF99CC", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Salinity", x ="Layer")

names(MP_abund)
rownames(MP_abund)
x<-tidyr::gather(MP_abund, variables, values, MP.RGC.v4_001193594:MP.RGC.v4_067950009)
x$values
x$variables
head(x)
names(x)
##Oxigen
class(x$oxygen)
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>% #excloure missing values
  select(c("sampling_station", "depth", "layer", "oxygen", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(oxygen), size = values)) +
  #geom_point(color = "#003399", alpha=0.7) +
  geom_point(color = "#3399FF", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Oxygen", x ="Layer")


##Fluorescencia
class(x$fluorescence)
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "fluorescence", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(fluorescence), size = values)) +
  geom_point(color = "#33FF99", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Fluorescence", x ="Layer")

##Turbidity
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "turbidity", "variables", "values")) %>%
  ggplot(aes(x=layer, y =turbidity, size = values)) +
  geom_point(color = "#006600", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Turbidity", x ="Layer")

##Nitrat
class(x$NO3)
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "NO3", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(NO3), size = values)) +
  geom_point(color = "#FF9966", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Nitrate concentration (µmol/Kg)", x ="Layer")


##Fosfat
#class(x$PO4)
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "PO4", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(PO4), size = values)) +
  geom_point(color = "#9966FF", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Phosphate concentration (µmol/Kg)", x ="Layer")


##Silici
#class(x$SIO2)
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "SIO2", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(SIO2), size = values)) +
  geom_point(color = "#Cc0066", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="SiO4 concentration (µmol/Kg)", x ="Layer")


##Latitud
x$layer<-gsub("epipelagic", "Epi", x$layer)
x$layer<-gsub("mesopelagic", "Meso", x$layer)
x$layer<-gsub("bathypelagic", "Bathy", x$layer)
x$layer<-factor(x$layer, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "latitude", "variables", "values")) %>%
  ggplot(aes(x=layer, y =latitude, size = values)) +
  geom_point(color = "#660066", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Latitude", x ="Layer")


names(MP_abund)
rownames(MP_abund)
x<-tidyr::gather(MP_abund, variables, values, MP.RGC.v4_001193594:MP.RGC.v4_067950009)
x$values
x$variables
head(x)
names(x)
x$layer<-gsub("epipelagic", "Epipelagic", x$layer)
x$layer<-gsub("mesopelagic", "Mesopelagic", x$layer)
x$layer<-gsub("bathypelagic", "Bathypelagic", x$layer)
x$layer<-factor(x$layer, levels=c("Epipelagic", "DCM", "Mesopelagic", "Bathypelagic" ))
class(x$values)
##Abundància mitjana
#vertical_palette <- c("epipelagic" = "#FFFF33", "DCM"= "#b7db70", "mesopelagic" = "#cce5ff", "bathypelagic"= "#cccccc")
vertical_palette <-c("#FFFF33","#b7db70","#cce5ff","#cccccc")
x %>% 
  data.frame()%>%
  #subset(sampling_station=="MH_019") %>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "longitude", "variables", "values")) %>%
  ggplot(aes(x=layer, y =as.numeric(values), fill=layer))+
  geom_boxplot(aes(fill=layer), size=0.2, alpha=0.8, notch = T) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=vertical_palette)+
  scale_y_log10()+ 
  theme_bw()+
  theme(legend.position="none")+
  labs(x="Layers", y ="LuxA gene mean abundance")


x$SumLuxA<-rowSums(MP_rgc_VP_t)
#x$SumLuxA<-Abundtotal_MP
##Matriu correlació
FQ<-as.data.frame((x[,c(12,7,8,15,14,18,17,21,28:30,33)]))
library(corrplot)
str(FQ)
FQ<-transform(FQ, fluorescence=as.numeric(fluorescence))
FQ<-transform(FQ, oxygen_concentration=as.numeric(oxygen_concentration))
FQ<-transform(FQ, values=as.numeric(values))
FQ<-transform(FQ, NO3=as.numeric(NO3))
FQ<-transform(FQ, PO4=as.numeric(PO4))
FQ<-transform(FQ, SIO2=as.numeric(SIO2))
FQ$rich<-vegan::specnumber(MP_rgc_VP_t) #afegir columna amb calcul de riquesa
FQ$div<-vegan::diversity(MP_rgc_VP_t) #afegir columna amb calcul de diversitat

colnames(FQ)<-c("Depth","Lat","Lon","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO2","LuxA","Rich","Div")

FQs<-as.data.frame(scale(FQ))
colnames(FQs)<-c("Depth","Lat","Lon","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO4","LuxA","Rich","Div")
FQcor<-cor(FQs, use="complete.obs")
corrplot.mixed(FQcor,lower="number",upper="ellipse",tl.col="black",tl.cex=0.8, number.cex=0.8)

##Matriu correlació estacio 19
library(corrplot)
FQ<-as.data.frame((x[,c(4,12,15,14,18,17,21,28:30,32)]))
str(FQ)
FQ$sampling_station <-gsub("[A-Za-z_]","", FQ$sampling_station) 
FQ<-transform(FQ, sampling_station=as.numeric(sampling_station))
FQ<-transform(FQ, fluorescence=as.numeric(fluorescence))
FQ<-transform(FQ, oxygen_concentration=as.numeric(oxygen_concentration))
FQ<-transform(FQ, values=as.numeric(values))
FQ<-transform(FQ, NO3=as.numeric(NO3))
FQ<-transform(FQ, PO4=as.numeric(PO4))
FQ<-transform(FQ, SIO2=as.numeric(SIO2))
#FQs<-as.data.frame(scale(FQ))
FQ$rich<-vegan::specnumber(MP_rgc_VP_t) #afegir columna amb calcul de riquesa
FQ$div<-vegan::diversity(MP_rgc_VP_t) #afegir columna amb calcul de diversitat

colnames(FQ)<-c("Station","Depth","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO2","LuxA","Rich","Div")
station_19_data<-subset(FQ, Station==19)
station_19_data<-na.omit(station_19_data[,c("Depth","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO2","LuxA","Rich","Div")])
cor_matrix_19<-cor(station_19_data[,c("Depth","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO2","LuxA","Rich","Div")])
par(mar=c(5,3,5,5))
corrplot(cor_matrix_19, method="circle", title= paste("Correlation Plot for st. 19"))
corrplot.mixed(cor_matrix_19,lower="number",upper="ellipse",tl.col="black",tl.cex=0.8, number.cex=0.8, title= paste("MH_019"))

##Matriu correlació estacions
library(corrplot)
FQ<-as.data.frame((x[,c(4,12,15,14,18,17,21,28:30,32)]))
str(FQ)
FQ$sampling_station <-gsub("[A-Za-z_]","", FQ$sampling_station) 
FQ<-transform(FQ, sampling_station=as.numeric(sampling_station))
FQ<-transform(FQ, fluorescence=as.numeric(fluorescence))
FQ<-transform(FQ, oxygen_concentration=as.numeric(oxygen_concentration))
FQ<-transform(FQ, values=as.numeric(values))
FQ<-transform(FQ, NO3=as.numeric(NO3))
FQ<-transform(FQ, PO4=as.numeric(PO4))
FQ<-transform(FQ, SIO2=as.numeric(SIO2))
#FQs<-as.data.frame(scale(FQ))
FQ$rich<-vegan::specnumber(MP_rgc_VP_t) #afegir columna amb calcul de riquesa
FQ$div<-vegan::diversity(MP_rgc_VP_t) #afegir columna amb calcul de diversitat

colnames(FQ)<-c("Station","Depth","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO2","LuxA","Rich","Div")
station_X_data<-subset(FQ, Station==141)
station_X_data<-na.omit(station_X_data[,c("Depth","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO2","LuxA","Rich","Div")])
cor_matrix_X<-cor(station_X_data[,c("Depth","Temp","Sal","Fl","O2","Turb","NO3","PO4","SiO2","LuxA","Rich","Div")])
#corrplot(cor_matrix_X, method="circle", title= paste("Correlation Plot for st.X"))
corrplot.mixed(cor_matrix_X,lower="number",upper="ellipse",tl.col="black",tl.cex=0.8, number.cex=0.8)

#Rectes regressió 
names(MP_abund)
rownames(MP_abund)
x<-tidyr::gather(MP_abund, variables, values, MP.RGC.v4_001193594:MP.RGC.v4_067950009)
x$values
x$variables
head(x)
names(x)
class(x$depth)
range(x$depth)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "depth", "layer", "variables", "values","bac_conc")) %>%
  #filter(bac_conc<1.5e+06)%>%
  ggplot(aes(y=as.numeric(bac_conc), x =depth))+
  geom_point()+
  geom_smooth(method="lm", se=F, color="#0000cc")+
  scale_y_log10()+
  scale_x_log10()+
  #scale_y_continuous(labels=function(x) format(x,scientific=TRUE))+ ##posar numeros en format cientific
  #scale_y_reverse()+
  theme_bw()+
  theme(legend.position="none")+
  labs(y="Bacteria cell abundance (cells/ml)", x ="Depth (m)")


#Calculate diversity and richness, and add them as variables to the env table
?vegan::specnumber
MP_rgc_VP_t<-MP_rgc_VP_t[,-1]
class(MP_rgc_VP_t)
str(MP_rgc_VP_t)
rownames(malaspina.codis)
colnames(malaspina.codis)
rownames(MP_rgc_VP_t)
colnames(MP_rgc_VP_t)
MP_rgc_VP_t<-MP_rgc_VP_t%>%
  as.data.frame()%>%
  mutate_if(is.character, as.numeric) ##posar chracter as numeric
view(MP_rgc_VP_t)
malaspina.codis$rich<-vegan::specnumber(MP_rgc_VP_t) #afegir columna amb calcul de riquesa
malaspina.codis$div<-vegan::diversity(MP_rgc_VP_t) #afegir columna amb calcul de diversitat
view(malaspina.codis)
rowSums(MP_rgc_VP_t)
malaspina.codis$sampling_station<-factor(malaspina.codis$sampling_station, levels=c("MH_019", "MH_030", "MH_044", "MH_049", "MH_063","MH_076", "MH_083","MH_092", "MH_101", "MH_120", "MH_141" ))

##Riquesa
ggplot(malaspina.codis, aes(x=depth,y=rich))+
  # geom_boxplot(aes(fill=Stations), size=0.2,alpha=0.8,outlier.size=0.3, notch = F)+
  theme_bw()+
  #geom_point()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_smooth(method="lm", se=F, color="#0000cc")+
  #scale_fill_manual(values=c(vertical_palette),name="")+
  facet_wrap(~sampling_station,ncol=6, scales="fixed")+
  theme(strip.background=element_rect(fill = '#f3f3f3'))+
  theme(panel.background = element_rect(fill = 'transparent'), plot.background = element_rect(fill='transparent', color=NA),legend.background = element_rect(fill='transparent'))+
  theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),axis.ticks=element_line(size=0.2),aspect.ratio=8/(7+sqrt(2)),text=element_text(size=12),axis.text.y=element_text(size=9),
        axis.title.y=element_text(margin=margin(0,10,0,0)),axis.title.x=element_text(margin=margin(5,0,0,0)),axis.text.x=element_text(size=9,hjust=1,vjust=0.5,angle=90),
        plot.title=element_text(size=12, face='bold', family="Times",margin=margin(0,0,20,0)))+
  labs(title="",x="Depth (m)",y="LuxA gene richness")+theme(legend.position="right")


##Diversitat
ggplot(malaspina.codis, aes(x=depth,y=div))+
  #geom_boxplot(aes(fill=values),size=0.2,alpha=0.8,outlier.size=0.3, notch = F)+
  theme_bw()+
  #geom_point()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_smooth(method="lm", se=F, color="red")+
  #scale_fill_manual(values=c(superpalette),name="")+
  facet_wrap(~sampling_station,ncol=6, scales="fixed")+
  theme(strip.background=element_rect(fill = '#f3f3f3'))+
  theme(panel.background = element_rect(fill = 'transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        axis.ticks=element_line(size=0.2),aspect.ratio=8/(7+sqrt(2)),
        text=element_text(size=12),axis.text.y=element_text(size=9),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x=element_text(margin=margin(5,0,0,0)),
        axis.text.x=element_text(size=9,hjust=1,vjust=0.5,angle=90),
        plot.title=element_text(size=12, face='bold', family="Times",
                                margin=margin(0,0,20,0)))+
  labs(title="",x="Depth (m)",y="LuxA gene diveristy")+theme(legend.position="none")

##Diversitat/Riquesa
ggplot(malaspina.codis, aes(x=depth,y=div/rich))+
  theme_bw()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_smooth(method="lm", se=F, color="green")+
  facet_wrap(~sampling_station,ncol=6, scales="fixed")+
  theme(strip.background=element_rect(fill = '#f3f3f3'))+
  theme(panel.background = element_rect(fill = 'transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        axis.ticks=element_line(size=0.2),aspect.ratio=8/(7+sqrt(2)),
        text=element_text(size=12),axis.text.y=element_text(size=9),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x=element_text(margin=margin(5,0,0,0)),
        axis.text.x=element_text(size=9,hjust=1,vjust=0.5,angle=90),
        plot.title=element_text(size=12, face='bold', family="Times",
                                margin=margin(0,0,20,0)))+
  labs(title="",x="Depth (m)",y="Diveristy/Richness")+theme(legend.position="none")



######Step1 1: Install and load required packages
install.packages("leaflet")
install.packages("leaflet.extras") # For additional features (optional)
# Load the required packages
library(leaflet)
library(leaflet.extras) # If you installed it for additional features

######Step 2: Prepare data

#MP_rgc_VP_t<-MP_rgc_VP_t[,-1]
MP_rgc_sum<-(MP_rgc_VP_t)
MP_rgc_sum<-MP_rgc_sum%>%
  as.data.frame()%>%
  mutate_if(is.character, as.numeric) ##posar chracter as numeric
MP_sum<-as.data.frame(rowSums(MP_rgc_sum))
MP <- cbind(MP, MP_sum)

library("writexl")
m.mapa <- read_excel("C:/Users/Desktop/TFM/malaspina/taula_mapa_MP.xlsx")
m.mapa<-as.data.frame((m.mapa[,c(1,2,3,4)]))
str(m.mapa)
print(unique(m.mapa$sampling_station))


######Step 3: Create the leaflet map with bubble points
# Create the leaflet map
world_map <- leaflet() %>%
  addTiles() # Adds the base world map layer

min_radius<-5
max_radius<-20

scaled_radius<-((m.mapa$sum_st-min(m.mapa$sum_st))/(max(m.mapa$sum_st)-min(m.mapa$sum_st)))*(max_radius-min_radius)+min_radius

world_map <- world_map %>%
  addTiles(tileOptions(provider="OpenStreetMap")) #%>%
addCircleMarkers(
  data = m.mapa,
  lat = ~lat_av,      #si la columna amb la latitud es diu Lat posar lat = ~Lat
  lng = ~long_av,     # si columna amb longitud es diu Lon, long = ~Lon
  #radius = ~sum_st,  #si vols que la mida sigui segons l'abundància si lo columna es diu abund radius= ~abund
  radius = scaled_radius,
  color = "blue",      # Bubble outline color
  fillColor = "red",   # Bubble fill color
  fillOpacity = 0.7,   # Bubble fill opacity
  stroke = TRUE,       # Show bubble outline
  weight = 1,
  label = ~as.character(sampling_station))#%>%         # Outline weight


legend_html <- sprintf(
  '<div class="legend-container" style="position: absolute; bottom: 10px; right: 10px; width:100px; border: 1px solid #ccc; background-color: #fff; padding: 5px;">
    <div class="legend-title" style="font-size: 12px;">Size Legend</div>
    <div class="legend-item">
      <div class="legend-circle" style="width: %dpx; height: %dpx; background-color: blue; border-radius: 50%%;"></div>
      <div class="legend-label" style="font-size: 10px;">%g - %s</div>
    </div>
    <div class="legend-item">
      <div class="legend-circle" style="width: %dpx; height: %dpx; background-color: blue; border-radius: 50%%;"></div>
      <div class="legend-label" style="font-size: 10px;">%g - %s</div>
    </div>
  </div>',
  min_radius, min_radius, min(m.mapa$sum_st), m.mapa$sampling_station[which.min(m.mapa$sum_st)],
  max_radius, max_radius, max(m.mapa$sum_st), m.mapa$sampling_station[which.max(m.mapa$sum_st)]
)

# Add the legend to the map
world_map <- world_map %>%
  addControl(
    html = legend_html,
    position = "bottomright",
    className = "custom-legend middle"
  )
#addCSSRules(".custom-legend{top:40%}")


# Print the map
world_map


########################## HOTMIX ################################
library(readxl)
HM_codis <- read.delim("C:/Users/Desktop/TFM/hotmix/Hotmix_metadata_metagenomas.txt")
HM_rgc_VP <- read_excel("C:/Users/Desktop/TFM/hotmix/Counts_Hotmix_exclo.xlsx") #num grans->NA
HM_fusion <- read_excel("C:/Users/Desktop/TFM/hotmix/Hotmix_fusion.xlsx")

HM_rgc_VP_t <- t(HM_rgc_VP) #transposar taula pq quedi com hotmix_codis
colnames(HM_rgc_VP_t) <- HM_rgc_VP_t[1,] #posar primera fila com a nom de columna
HM_rgc_VP_t<-HM_rgc_VP_t[-1,] #eliminar primera fila pq ja esta com a nom de columnes

HM_rgc_VP_t<-data.frame(HM_rgc_VP_t)

names(HM_fusion)
rownames(HM_fusion)


#facet_wrap(~Stations)
HM_fusion<-data.frame(HM_fusion)
x<-tidyr::gather(HM_fusion, variables, values, HM_MG11_290740:HM_MG9_87741)
x$values
x$variables
head(x)
names(x)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
superpalette<-c("#8D4C6A", #cyano-cyano
                "#c2abc9", #gamma-unk
                "#b27300", #bacte-unk
                "#FFDF26", #unk-unk (0.2 no)
                "#626591", #alpha-cladeI
                "#0000ff", #alpha-cladeII (0.2)
                "#d8e6f3", #alpha-cladeIII
                "#66b2b2", #alpha-rhodo
                "#00004c", #alpha-stappia
                "#005b96", #alpha-AEGEAN (0.2)
                "#006666", #alpha-hypho (0.2 no)
                "#377EB8", #alpha-parv (0.2 no)
                "#00ff99", #alpha-sphingo
                "#4DAF4A", #gamma-comamo
                "#74ee15", #gamma-methyl
                "#028900", #gamma-burkhol
                "#740001", #gamma-oxalo  (0.2 no)
                "#ff2032", #gamma-morax
                "#4F507F", #gamma-pseudomona
                "#99cc99", #gamma-litorico
                "#6073b6", #gamma-spongiibact  (0.2 no)
                "#fe9360", #bacter-flavobacteria
                "#ff00b4", #bacter-NS9
                "#f3fc64", #Bacter-Crocinito
                "#b12ef7", #bacter-cryomo
                "#a7dbd8", #bacter-saprospir  (0.2 no)
                "#ea80fc", #bacter-chitinoph
                "#310e59", #bacter-NS11-12
                "#B15A7A", #bacter-cyclobact
                "#ff2e65", #bacter-spirosom (0.2)
                "#e6e8ea", #bacter-hymeno (0.2 no)
                "#FF7F00", #rhodothermia-balneolacae
                "#dab600", #actino-sporichthya
                "#ff80ed", #actino-microbact
                "#F15E75", #verruco-rubrital (0.2 no)
                "#FFFF33", #verruco-punicei
                "#ff0000", #amatimona-amatimona
                "#4169E1", #bacillii-planococc (0.2 no)
                "#0000CD") #z_other
palette<-colorRampPalette(cbbPalette)
palette_2<-colorRampPalette(superpalette)

##Abundàncies totals
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>%
  data.frame()%>%
  na.exclude()%>% #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "variables", "values")) %>%
  #filter(values<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(values), fill = variables)) +
  #geom_point() + geom_line()+
  geom_bar(position ="stack", stat="identity")+#graph amb abund totals (position ="stack")
  scale_fill_manual(values=palette_2(90))+    
  #scale_y_log10() +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(x="Layers", y ="LuxA gene total abundance")


##Abundàncies relatives
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "variables", "values")) %>%
  #filter(as.numeric(values)<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(values), fill = variables)) +
  geom_bar(position ="fill", stat="identity")+ #graph amb abund relatives (position="fill") 
  scale_fill_manual(values=palette_2(90))+
  # scale_y_log10() +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(x="Layers", y ="LuxA gene relative abundance")


##Gràfic depth/long 
names(HM_fusion)
rownames(HM_fusion)
#facet_wrap(~Stations)
HM_fusion<-data.frame(HM_fusion)
x<-tidyr::gather(HM_fusion, variables, values, HM_MG11_290740:HM_MG9_87741)
x$values
x$variables
head(x)
names(x)

x$LAYER<-factor(x$LAYER, levels=c("Epipelagic", "DCM", "Mesopelagic", "Bathypelagic" ))

##Distribució vertical
vertical_palette <- c("#FFFF33","#b7db70","#cce5ff","#cccccc")
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "Lon", "variables", "values")) %>%
  #filter(values<2)%>%
  ggplot(aes(x=as.numeric(Lon), y =DEPTH, size = as.numeric(values), color = as.character(Stations), shape = LAYER)) +
  geom_rect(inherit.aes=F, aes(xmin=-Inf, xmax=Inf, ymin=0, max = 200), fill = "#ffff99", alpha =0.02)+
  geom_rect(inherit.aes=F, aes(xmin=-Inf, xmax=Inf, ymin=200, max = 1000), fill = "#cce5ff", alpha =0.02)+
  geom_rect(inherit.aes=F, aes(xmin=-Inf, xmax=Inf, ymin=1000, max = 2500), fill = "#cccccc", alpha =0.02)+
  geom_point(alpha=0.7) +
  scale_shape_manual(values = c("Epipelagic"=16,"DCM"= 13, "Mesopelagic"= 16, "Bathypelagic"=16))+
  scale_size_continuous(range= c(1,10),name=NULL)+
  scale_y_reverse()+
  #scale_y_continuous(trans=trans_reverser('log10'))+
  scale_color_manual(values=palette_2(11))+
  theme_bw()+
  theme(legend.position="right")+
  labs(x="Longitude", y ="Depth (m)", size= "Size legend", color="Stations", shape="Layer")
#guides(LAYER=guide_legend (title=NULL),as.character(Stations)=guide_legend (title=NULL))

##Temperatura
#graph parameters
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "T090C", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =T090C, size = values)) +
  geom_point(color = "#FF6633", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Temperature (ºC)", x ="Layer")

##Salinitat
#graph parameters
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "Sal00", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =Sal00, size = values)) +
  geom_point(color = "#FF99CC", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Salinity (PSU)", x ="Layer")


##Oxigen
#graph parameters
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "O2_CTD_umol.Kg", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =O2_CTD_umol.Kg, size = values)) +
  geom_point(color = "#3399FF", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Oxygen concentration (µmol/Kg)", x ="Layer")


names(HM_fusion)
rownames(HM_fusion)
#facet_wrap(~Stations)
HM_fusion<-data.frame(HM_fusion)
x<-tidyr::gather(HM_fusion, variables, values, HM_MG11_290740:HM_MG9_87741)
x$values
x$variables
head(x)
names(x)

##Fluorescencia
class(x$FlSPSbe)
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("sampling_station", "DEPTH", "LAYER", "FlSPSbe", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =FlSPSbe, size = values)) +
  geom_point(color = "#33FF99", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~sampling_station)+
  labs(y="Fluorescence", x ="Layer")

##Turbidity
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "SeaTurbMtr", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =SeaTurbMtr, size = values)) +
  geom_point(color = "#006600", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Turbidity (FTU)", x ="Layer")

##Nitrat
class(x$NO3_umol.kg)
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "NO3_umol.kg", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(NO3_umol.kg), size = values)) +
  geom_point(color = "#FF9966", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Nitrate concentration (µmol/Kg)", x ="Layer")


##Fosfat
class(x$PO4_umol.kg)
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "PO4_umol.kg", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(PO4_umol.kg), size = values)) +
  geom_point(color = "#9966FF", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Phosphate concentration (µmol/Kg)", x ="Layer")


##Silici
class(x$SIO2_umol.kg)
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "SIO2_umol.kg", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(SIO2_umol.kg), size = values)) +
  geom_point(color = "#Cc0066", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Silica concentration (µmol/Kg)", x ="Layer")

##Bak
class(x$Bacteria_.cells.ml.)
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "Bacteria_.cells.ml.", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(Bacteria_.cells.ml.), size = values)) +
  geom_point(color = "#cc9900", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="Bacteria concentration (cells/ml)", x ="Layer")


##TOC
class(x$TOC_.umolC.L.)
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
names(x)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "TOC_.umolC.L.", "variables", "values")) %>%
  filter(values<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(TOC_.umolC.L.), size = values)) +
  geom_point(color = "#66cc99", alpha=0.7) +
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Stations)+
  labs(y="TOC concentration (µmol C/L)", x ="Layer")


names(HM_fusion)
rownames(HM_fusion)
#facet_wrap(~Stations)
HM_fusion<-data.frame(HM_fusion)
x<-tidyr::gather(HM_fusion, variables, values, HM_MG11_290740:HM_MG9_87741)
x$values
x$variables
head(x)
names(x)
x$LAYER<-gsub("Epipelagic", "Epi", x$LAYER)
x$LAYER<-gsub("Mesopelagic", "Meso", x$LAYER)
x$LAYER<-gsub("Bathypelagic", "Bathy", x$LAYER)
x$LAYER<-factor(x$LAYER, levels=c("Epi", "DCM", "Meso", "Bathy" ))
class(x$values)
##Abundància mitjana
#vertical_palette <- c("epipelagic" = "#FFFF33", "DCM"= "#b7db70", "mesopelagic" = "#cce5ff", "bathypelagic"= "#cccccc")
vertical_palette <-c("#FFFF33","#b7db70","#cce5ff","#cccccc")
x %>% 
  data.frame()%>%
  #subset(sampling_station=="MH_019") %>%
  na.exclude()%>%  #excloure missing values
  select(c("Stations", "DEPTH", "LAYER", "Lon", "variables", "values")) %>%
  #filter(values<2)%>%
  ggplot(aes(x=LAYER, y =as.numeric(values), fill=LAYER))+
  #geom_boxplot(color="black", fill="orange", alpha=0.2)+
  geom_boxplot(aes(fill=LAYER), size=0.2, alpha=0.8, notch = T) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=vertical_palette)+
  scale_y_log10()+ 
  theme_bw()+
  theme(legend.position="none")+
  labs(x="Layers", y ="LuxA gene mean abundance")

##Matriu correlació
FQ<-as.data.frame((x[,c(8,11:15,20,29:31,48,57,60:62,65)]))#%>%
# filter(values<2)
#library(corrplot)
str(FQ)
FQ<-transform(FQ, NO3_umol.kg=as.numeric(NO3_umol.kg))
FQ<-transform(FQ, PO4_umol.kg=as.numeric(PO4_umol.kg))
FQ<-transform(FQ, SIO2_umol.kg=as.numeric(SIO2_umol.kg))
FQ<-transform(FQ, values=as.numeric(values))
FQ$rich<-vegan::specnumber(HM_rgc_VP_t) #afegir columna amb calcul de riquesa
FQ$div<-vegan::diversity(HM_rgc_VP_t)
FQs<-as.data.frame(scale(FQ))
colnames(FQs)<-c("Depth","Temp","Sal","Fl","o2","Turb","pH","NO3","PO4","SiO2","TOC","Bac","N","C","C/N","LuxA","Rich","Div")
FQcor<-cor(FQs, use="complete.obs")
corrplot.mixed(FQcor,lower="number",upper="ellipse",tl.col="black",tl.cex=0.8, number.cex=0.8)


#Rectes regressió 
names(HM_fusion)
rownames(HM_fusion)
#facet_wrap(~Stations)
HM_fusion<-data.frame(HM_fusion)
x<-tidyr::gather(HM_fusion, variables, values, HM_MG11_290740:HM_MG9_87741)
head(x)
names(x)
class(x$values)
x %>% 
  data.frame()%>%
  na.exclude()%>%  #excloure missing values
  select(c("DEPTH", "variables", "values", "Bacteria_.cells.ml.")) %>%
  # filter(values<2)%>%
  ggplot(aes(y=as.numeric(Bacteria_.cells.ml.), x=DEPTH))+
  geom_point()+
  geom_smooth(method="lm", se=F, color="#0000cc")+
  #scale_y_reverse()+
  scale_y_log10()+ 
  scale_x_log10()+
  theme_bw()+
  theme(legend.position="none")+
  labs(y="Bacteria cell abundance (cells/ml)", x="Depth (m)")


#Calculate diversity and richness (number of ASVs), and add them as variables to the env table
?vegan::specnumber
class(HM_rgc_VP_t)
str(HM_rgc_VP_t)
HM_rgc_VP_t<-HM_rgc_VP_t%>%
  as.data.frame()%>%
  mutate_if(is.character, as.numeric)#%>% ##posar character as numeric
####HM_rgc_VP_t<-na.omit(HM_rgc_VP_t)
#filter(HM_MG11_290740<2)
#HM_rgc_VP_t_filtered<-HM_rgc_VP_t%>% ##per treure valors erronis
# filter_all(all_vars(.<2))
#HM_rgc_VP_t[HM_rgc_VP_t>2]<-NA #
HM_rgc_VP_t[is.na(HM_rgc_VP_t)]<-0
View(HM_rgc_VP_t)
HM_rgc_VP_t<-HM_rgc_VP_t[-1,]

HM_codis$rich<-specnumber(HM_rgc_VP_t)
HM_codis$rich<-vegan::specnumber(HM_rgc_VP_t) #afegir columna amb calcul de riquesa
HM_codis$div<-vegan::diversity(HM_rgc_VP_t) #afegir columna amb calcul de diversitat
view(HM_codis)
rownames(HM_codis)
colnames(HM_codis)
rowSums(HM_rgc_VP_t)
HM_codis$Stations<-factor(HM_codis$Stations, levels=c("St.2", "St.9", "St.17", "St.22", "St.25" ))
#Riquesa
ggplot(HM_codis, aes(x=depth,y=rich))+
  # geom_boxplot(aes(fill=Stations), size=0.2,alpha=0.8,outlier.size=0.3, notch = F)+
  theme_bw()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_smooth(method="lm", se=F, color="#0000cc")+
  facet_wrap(~Stations,ncol=5, scales="fixed")+
  theme(strip.background=element_rect(fill = '#f3f3f3'))+
  theme(panel.background = element_rect(fill = 'transparent'), plot.background = element_rect(fill='transparent', color=NA),legend.background = element_rect(fill='transparent'))+
  theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),axis.ticks=element_line(size=0.2),aspect.ratio=8/(7+sqrt(2)),text=element_text(size=12),axis.text.y=element_text(size=9),
        axis.title.y=element_text(margin=margin(0,10,0,0)),axis.title.x=element_text(margin=margin(5,0,0,0)),axis.text.x=element_text(size=9,hjust=1,vjust=0.5,angle=90),
        plot.title=element_text(size=12, face='bold', family="Times",margin=margin(0,0,20,0)))+
  labs(title="",x="Depth (m)",y="luxA gene richness")+theme(legend.position="right")

#Diversitat
ggplot(HM_codis, aes(x=depth,y=div))+
  theme_bw()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_smooth(method="lm", se=F, color="red")+
  facet_wrap(~Stations,ncol=5, scales="fixed")+
  theme(strip.background=element_rect(fill = '#f3f3f3'))+
  theme(panel.background = element_rect(fill = 'transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        axis.ticks=element_line(size=0.2),aspect.ratio=8/(7+sqrt(2)),
        text=element_text(size=12),axis.text.y=element_text(size=9),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x=element_text(margin=margin(5,0,0,0)),
        axis.text.x=element_text(size=9,hjust=1,vjust=0.5,angle=90),
        plot.title=element_text(size=12, face='bold', family="Times",
                                margin=margin(0,0,20,0)))+
  labs(title="",x="Depth (m)",y="luxA gene diveristy")+theme(legend.position="none")


##Diversitat/Riquesa
ggplot(HM_codis, aes(x=depth,y=div/rich))+
  theme_bw()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_smooth(method="lm", se=F, color="green")+
  facet_wrap(~Stations,ncol=5, scales="fixed")+
  theme(strip.background=element_rect(fill = '#f3f3f3'))+
  theme(panel.background = element_rect(fill = 'transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        axis.ticks=element_line(size=0.2),aspect.ratio=8/(7+sqrt(2)),
        text=element_text(size=12),axis.text.y=element_text(size=9),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x=element_text(margin=margin(5,0,0,0)),
        axis.text.x=element_text(size=9,hjust=1,vjust=0.5,angle=90),
        plot.title=element_text(size=12, face='bold', family="Times",
                                margin=margin(0,0,20,0)))+
  labs(title="",x="Depth (m)",y="Diveristy/Richness")+theme(legend.position="none")


#Mapa bubbleplot
#Step 1: Load the required packages
library(leaflet)
library(leaflet.extras) # If you installed it for additional features

######Step 2: Prepare data
#(taula amb les coordenades de cada estació i una columna amb l'abundància a cada estació)
HM_rgc_sum<-(HM_rgc_VP_t)
HM_rgc_sum<-HM_rgc_sum%>%
  as.data.frame()%>%
  mutate_if(is.character, as.numeric) ##posar chracter as numeric
HM_sum<-as.data.frame(rowSums(HM_rgc_sum))

library("writexl")
write_xlsx(HM_sum, "C:/Users/Desktop/TFM/hm_table.xlsx")
hm.mapa <- read_excel("C:/Users/Desktop/TFM/hotmix/taula_mapa_hm.xlsx")

hm.mapa<-as.data.frame((hm.mapa[,c(1,2,3,4)]))
str(hm.mapa)
print(unique(hm.mapa$Stations))

######Step 3: Create the leaflet map with bubble points
# Create the leaflet map
world_map <- leaflet() %>%
  addTiles() # Adds the base world map layer

min_radius<-5
max_radius<-20

scaled_radius<-((hm.mapa$rowSums-min(hm.mapa$rowSums))/(max(hm.mapa$rowSums)-min(hm.mapa$rowSums)))*(max_radius-min_radius)+min_radius

world_map <- world_map %>%
  addTiles() %>%
  addCircleMarkers(
    data = hm.mapa,
    lat = ~Lat,      
    lng = ~Lon,     
    #radius = ~sum_st,  #queden massa petits
    radius = scaled_radius,
    color = "blue",      # Bubble outline color
    fillColor = "darkgreen",   # Bubble fill color
    fillOpacity = 0.7,   # Bubble fill opacity
    stroke = TRUE,       # Show bubble outline
    weight = 1,
    label = ~as.character(Stations),
    labelOptions=labelOptions(noHide=TRUE, style= "font-size:10px"))#%>%         # Outline weight


legend_html <- sprintf(
  '<div class="legend-container" style="position: absolute; bottom: 10px; right: 10px; width:100px; border: 1px solid #ccc; background-color: #fff; padding: 5px;">
    <div class="legend-title" style="font-size: 12px;">Size Legend</div>
    <div class="legend-item">
      <div class="legend-circle" style="width: %dpx; height: %dpx; background-color: green; border-radius: 50%%;"></div>
      <div class="legend-label" style="font-size: 10px;">%g - %s</div>
    </div>
    <div class="legend-item">
      <div class="legend-circle" style="width: %dpx; height: %dpx; background-color: green; border-radius: 50%%;"></div>
      <div class="legend-label" style="font-size: 10px;">%g - %s</div>
    </div>
  </div>',
  min_radius, min_radius, min(hm.mapa$rowSums), hm.mapa$Stations[which.min(hm.mapa$rowSums)],
  max_radius, max_radius, max(hm.mapa$rowSums), hm.mapa$Stations[which.max(hm.mapa$rowSums)]
)

# Add the legend to the map
world_map <- world_map %>%
  addControl(
    html = legend_html,
    position = "bottomright",
    className = "custom-legend middle"
  )

# Print the map
world_map
