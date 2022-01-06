
##########################################################################
######                                                              ######
######   	            stochastic character mapping                  ######
######                                                              ######
##########################################################################

library(phytools)
library(ape)
library(ggplot2)
library(cowplot)
library(coda)
library(gridExtra)
library(ggpubr)
library(tidyverse)

## Full Dataset

## Read in  full trait datasets
Full.Strict.Dataset <- read.csv("FullStrict_TraitAnalysis.csv",row.names=1)
Full.Relaxed.Dataset <- read.csv("FullRelaxed_TraitAnalysis.csv",row.names=1)

FullTree <- read.nexus(file = "FullTree.tre")


#Read in ecological and morphological (toepad) trait dataset (Supplementary Table 1) 
#Name rows with concatenated genus_species column

#Subset some characters for subsequent analyses

Full.Strict.Arboreal <- as.factor(setNames(Full.Strict.Dataset[,1],rownames(Full.Strict.Dataset))) #Strict arboreal Full dataset
Full.Toepads <- as.factor(setNames(Full.Strict.Dataset[,2],rownames(Full.Strict.Dataset))) #toepads ful ldataset

#Assign prior probability on root node for non-arboreality

arbor.pi<-setNames(c(1,0),c("0","1"))
arbor.pi


#fit two Mk models

#fitARD.strict.arbor.full<-fitMk(FullTree,Full.Strict.Arboreal,pi=arbor.pi,model="ARD")
#AIC(fitARD.strict.arbor.full)
#
#fitER.strict.arbor.full<<-fitMk(FullTree,Full.Strict.Arboreal,pi=arbor.pi,model="ER")
#AIC(fitER.strict.arbor.full)
#
#Compare AIC values and determine best model

#aic.strict.arbor.full<-setNames(sapply(list(fitARD.strict.arbor.full,fitER.strict.arbor.full),AIC), c("fitARD","fitER"))
#aic.strict.arbor.full
#aic.w(aic.strict.arbor.full)

## Full Strict Arboreal reconstructions

Full.Strict.Arboreal.simmap.ARD<-make.simmap(FullTree,Full.Strict.Arboreal,nsim=500,pi=arbor.pi,model="ARD")
Full.Strict.Arboreal.simmap.ER<-make.simmap(FullTree,Full.Strict.Arboreal,nsim=500,pi=arbor.pi,model="ER")

#plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Full.Strict.Arboreal.simmap.ARD <-densityMap(Full.Strict.Arboreal.simmap.ARD,states=levels(Full.Strict.Arboreal)[1:2],plot=FALSE)
colors.1<-setMap(Obj.Full.Strict.Arboreal.simmap.ARD,c("gray84","lightskyblue"))
plot(colors.1,Obj.Full.Strict.Arboreal.simmap.ARD,outline=F,fsize=c(0.7,0.9),ftype="off",type = "fan",offset=1, lwd=1)
title(main = "Full Strict Arbor ARD 500")
ape::tiplabels(pch=21, tip=ii, cex=1.5, col='black',bg="grey")

## Padbearing clade arc labels
#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(FullTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(FullTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=findMRCA(FullTree,c("Typhlosaurus_braini","Carlia_dogare")),
                ln.offset=1.06,lab.offset=1.09, mark.node=FALSE)

## Families with arboreal species

##Agamidae
#arc.cladelabels(text="Agamidae",node=findMRCA(FullTree,c("Saara_asmussi","Pseudocalotes_brevipes")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "vertical", mark.node=FALSE)
##Gekkonidae
#arc.cladelabels(text="Gekkonidae",node=findMRCA(FullTree,c("Lepidodactylus_novaeguineae","Phelsuma_lineata")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Anguidae
#arc.cladelabels(text="Anguidae",node=findMRCA(FullTree,c("Abronia_campbelli","Celestus_enneagrammus")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "vertical", mark.node=FALSE)
##Chamaeleonidae
#arc.cladelabels(text="Chamaeleonidae",node=findMRCA(FullTree,c("Trioceros_rudis","Palleon_nasus")),
#                ln.offset=1.05,lab.offset=1.07,orientation = "vertical", mark.node=FALSE)
##Corytophanidae
#arc.cladelabels(text="Corytophanidae",node=findMRCA(FullTree,c("Basiliscus_plumifrons","Laemanctus_longipes")),
#                ln.offset=1.09,lab.offset=1.12, orientation = "vertical",mark.node=FALSE)
##Dactyloidae
#arc.cladelabels(text="Dactyloidae",node=findMRCA(FullTree,c("Anolis_bonairensis","Anolis_poecilopus")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Diplodactylidae
#arc.cladelabels(text="Diplodactylidae",node=findMRCA(FullTree,c("Pseudothecadactylus_lindneri","Diplodactylus_pulcher")),
#                ln.offset=1.03,lab.offset=1.06,mark.node=FALSE)
##Eublepharidae
#arc.cladelabels(text="Eublepharidae",node=findMRCA(FullTree,c("Aeluroscalabotes_felinus","Hemitheconyx_caudicinctus")),
#                ln.offset=1.15,lab.offset=1.18,orientation = "horizontal", mark.node=FALSE)
##Gekkonidae
#arc.cladelabels(text="Gekkonidae",node=findMRCA(FullTree,c("Lepidodactylus_novaeguineae","Phelsuma_lineata")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##2647 
##getDescendants(FullTree, 2647, curr=NULL)
#
##Gerrhosauridae
#arc.cladelabels(text="Gerrhosauridae",node=findMRCA(FullTree,c("Broadleysaurus_major","Zonosaurus_bemaraha")),
#                ln.offset=1.03,lab.offset=1.05,orientation = "vertical", mark.node=FALSE)
##Iguanidae
#arc.cladelabels(text="Iguanidae",node=findMRCA(FullTree,c("Dipsosaurus_dorsalis","Ctenosaura_palearis")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Lacertidae
#arc.cladelabels(text="Lacertidae",node=findMRCA(FullTree,c("Psammodromus_hispanicus","Darevskia_caucasica")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Leiosauridae
#arc.cladelabels(text="Leiosauridae",node=findMRCA(FullTree,c("Enyalius_bilineatus","Leiosaurus_paronae")),
#                ln.offset=1.21,lab.offset=1.24, orientation = "horizontal",mark.node=FALSE)
#
##Liolaemidae
#arc.cladelabels(text="Liolaemidae",node=findMRCA(FullTree,c("Ctenoblepharys_adspersa","Liolaemus_ornatus")),
#                ln.offset=1.03,lab.offset=1.06, orientation = "horizontal",mark.node=FALSE)
#
##Opluridae
#arc.cladelabels(text="Opluridae",node=findMRCA(FullTree,c("Chalarodon_madagascariensis","Oplurus_grandidieri")),
#                ln.offset=1.15,lab.offset=1.18,orientation = "horizontal", mark.node=FALSE)
##Phrynosomatidae
#arc.cladelabels(text="Phrynosomatidae",node=findMRCA(FullTree,c("Callisaurus_draconoides","Sceloporus_stejnegeri")),
#                ln.offset=1.09,lab.offset=1.12,orientation = "horizontal", mark.node=FALSE)
##Phyllodactylidae
#arc.cladelabels(text="Phyllodactylidae",node=findMRCA(FullTree,c("Thecadactylus_rapicauda","Tarentola_gigas")),
#                ln.offset=1.05,lab.offset=1.08,orientation = "horizontal", mark.node=FALSE)
##Polychrotidae
#arc.cladelabels(text="Polychrotidae",node=findMRCA(FullTree,c("Polychrus_gutturosus","Polychrus_marmoratus")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Scincidae
#arc.cladelabels(text="Scincidae",node=3052, ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Sphaerodactylidae
#arc.cladelabels(text="Sphaerodactylidae",node=findMRCA(FullTree,c("Sphaerodactylus_ocoae","Quedenfeldtia_moerens")),
#                ln.offset=1.03,lab.offset=1.06, orientation = "horizontal", mark.node=FALSE)
##Teiidae
#arc.cladelabels(text="Teiidae",node=findMRCA(FullTree,c("Pholidoscelis_fuscatus","Callopistes_maculatus")),
#                ln.offset=1.05,lab.offset=1.08, mark.node=FALSE)
##Tropiduridae
#arc.cladelabels(text="Tropiduridae",node=findMRCA(FullTree,c("Stenocercus_ochoai","Tropidurus_cocorobensis")),
#                ln.offset=1.05,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Varanidae
#arc.cladelabels(text="Varanidae",node=findMRCA(FullTree,c("Varanus_griseus","Varanus_acanthurus")),
#                ln.offset=1.09,lab.offset=1.12, mark.node=FALSE)
#


Obj.Full.Strict.Arboreal.simmap.ER <-densityMap(Full.Strict.Arboreal.simmap.ER,states=levels(Full.Strict.Arboreal)[1:2],plot=FALSE)
colors.2<-setMap(Obj.Full.Strict.Arboreal.simmap.ER,c("gray84","lightskyblue"))
plot(colors.2,Obj.Full.Strict.Arboreal.simmap.ER,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Full Strict Arbor ER 500")

#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(FullTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(FullTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=findMRCA(FullTree,c("Typhlosaurus_braini","Carlia_dogare")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Label padbearing species on tips
Full.toepad.tips<-read.csv("FullStrict_TraitAnalysis.csv")
Full.toepad.tips <- filter(Full.toepad.tips, toepads == "1")
Full.Toepad.Species <- Full.toepad.tips$Concatenated
ii<-sapply(Full.Toepad.Species,grep,FullTree$tip.label)
ii
ape::tiplabels(pch=21, tip=ii, cex=1.5, col='black',bg="grey")

#Examine the distribution of estimated origins and losses of arboreality in the strict arboreality dataset

Full.Strict.Arboreal.ARD.HPD<-describe.simmap(Full.Strict.Arboreal.simmap.ARD)
Full.Strict.Arboreal.ER.HPD<-describe.simmap(Full.Strict.Arboreal.simmap.ER)

#Summary statistics
Full.Strict.Arboreal.simmap.ARD.density <- density(Full.Strict.Arboreal.simmap.ARD) #Summary statistics
Full.Strict.Arboreal.simmap.ER.density <- density(Full.Strict.Arboreal.simmap.ER) #Summary statistics

#Strict arboreality frequency distribution
Full.Strict.Arbor.ARD.origins.losses<- Full.Strict.Arboreal.ARD.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Full.Strict.Arbor.ARD.origins.losses) <- c("Gains", "Losses") #rename columns
Full.Strict.Arbor.ARD.origins.losses.df <- as.data.frame(Full.Strict.Arbor.ARD.origins.losses)
#gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Full.Strict.Arbor.ARD <- ggplot(Full.Strict.Arbor.ARD.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_classic()+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+ggtitle("Gains.Distribution.Full.Strict.Arbor.ARD")+
  geom_text(x = 12,  y = 10, 
            label = "", 
            colour = "black")+
  geom_segment(x = 102, xend = 131, 
               y = 31, yend = 31,
               colour = "black")+scale_y_continuous(limits = c(0,33), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320)) 
#losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Full.Strict.Arbor.ARD <- ggplot(Full.Strict.Arbor.ARD.origins.losses.df, aes(x=Losses))+ggtitle("Losses.Distribution.Full.Strict.Arbor.ARD")+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_classic()+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 198, xend = 261, 
               y = 19, yend = 19,
               colour = "black")+scale_y_continuous(limits = c(0,20.5), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320)) 

#Put both histograms together
Full.Strict.Arbor.ARD.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Strict.Arbor.ARD, Losses.Distribution.Full.Strict.Arbor.ARD,labels = c("A", "B", ncol = 1, nrow = 1))


#Strict arboreality frequency distribution
Full.Strict.Arbor.ER.origins.losses<- Full.Strict.Arboreal.ER.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Full.Strict.Arbor.ER.origins.losses) <- c("Gains", "Losses") #rename columns
Full.Strict.Arbor.ER.origins.losses.df <- as.data.frame(Full.Strict.Arbor.ER.origins.losses)
#gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Full.Strict.Arbor.ER <- ggplot(Full.Strict.Arbor.ER.origins.losses.df, aes(x=Gains))+ggtitle("Gains.Distribution.Full.Strict.Arbor.ER")+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_classic()+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 152, xend = 179, 
               y = 38, yend = 38,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320)) 
#losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Full.Strict.Arbor.ER <- ggplot(Full.Strict.Arbor.ER.origins.losses.df, aes(x=Losses))+ggtitle("Losses.Distribution.Full.Strict.Arbor.ER")+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_classic()+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 66, xend = 94, 
               y = 33, yend = 33,
               colour = "black")+scale_y_continuous(limits = c(0,35), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320)) 

#Put both histograms together
Full.Strict.Arbor.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Strict.Arbor.ER, Losses.Distribution.Full.Strict.Arbor.ER,labels = c("A", "B", ncol = 1, nrow = 1))

# ARD and ER together
Full.Strict.Arbor.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Strict.Arbor.ARD, Losses.Distribution.Full.Strict.Arbor.ARD, Gains.Distribution.Full.Strict.Arbor.ER, Losses.Distribution.Full.Strict.Arbor.ER,labels = c("A", "B", "C", "D", ncol = 2, nrow = 2))


### all
#Full.Strict.Arbor.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Strict.Arbor.ARD, Losses.Distribution.Full.Strict.Arbor.ARD,Gains.Distribution.Full.Strict.Arbor.ER, Losses.Distribution.Full.Strict.Arbor.ER,labels = c("A", "B", "C", "D", ncol = 2, nrow = 2))

##########################################################################
######                                                              ######
######   	          Full Relaxed Arboreality                        ######
######                                                              ######
##########################################################################

Full.Relaxed.Dataset <- read.csv("FullRelaxed_TraitAnalysis.csv",row.names=1)

#Subset some characters for subsequent analyses

Full.Relaxed.Arboreal <- as.factor(setNames(Full.Relaxed.Dataset[,1],rownames(Full.Relaxed.Dataset))) #Relaxed arboreal Full dataset
Full.Toepads <- as.factor(setNames(Full.Relaxed.Dataset[,2],rownames(Full.Relaxed.Dataset))) #toepads ful ldataset

#Assign prior probability on root node for non-arboreality

#arbor.pi<-setNames(c(1,0),c("not.arboreal","arboreal"))

#fit two Mk models

#fitARD.Relaxed.arbor.full<-fitMk(FullTree,Full.Relaxed.Arboreal,pi=arbor.pi,model="ARD")
#AIC(fitARD.Relaxed.arbor.full)
#
#fitER.Relaxed.arbor.full<<-fitMk(FullTree,Full.Relaxed.Arboreal,pi=arbor.pi,model="ER")
#AIC(fitER.Relaxed.arbor.full)

#Compare AIC values and determine best model

#aic.Relaxed.arbor.full<-setNames(sapply(list(fitARD.Relaxed.arbor.full,fitER.Relaxed.arbor.full),AIC), c("fitARD","fitER"))
#aic.Relaxed.arbor.full
#aic.w(aic.Relaxed.arbor.full)

## Full Relaxed Arboreal reconstructions

Full.Relaxed.Arboreal.simmap.ARD<-make.simmap(FullTree,Full.Relaxed.Arboreal,nsim=500,pi=arbor.pi,model="ARD")
Full.Relaxed.Arboreal.simmap.ER<-make.simmap(FullTree,Full.Relaxed.Arboreal,nsim=500,pi=arbor.pi,model="ER")

#plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Full.Relaxed.Arboreal.simmap.ARD <-densityMap(Full.Relaxed.Arboreal.simmap.ARD,states=levels(Full.Relaxed.Arboreal)[1:2],plot=FALSE)
colors.3<-setMap(Obj.Full.Relaxed.Arboreal.simmap.ARD,c("gray84","lightskyblue"))
plot(colors.3,Obj.Full.Relaxed.Arboreal.simmap.ARD,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Full Relaxed Arbor ARD 500")

##Major padbearing clades

#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(FullTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(FullTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=findMRCA(FullTree,c("Typhlosaurus_braini","Carlia_dogare")),
                ln.offset=1.06,lab.offset=1.09, mark.node=FALSE)

#Label padbearing species on tips
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")


Obj.Full.Relaxed.Arboreal.simmap.ER <-densityMap(Full.Relaxed.Arboreal.simmap.ER,states=levels(Full.Relaxed.Arboreal)[1:2],plot=FALSE)
colors.4<-setMap(Obj.Full.Relaxed.Arboreal.simmap.ER,c("gray84","lightskyblue"))
plot(colors.4,Obj.Full.Relaxed.Arboreal.simmap.ER,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Full Relaxed Arbor ER 500")

#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(FullTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(FullTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=findMRCA(FullTree,c("Typhlosaurus_braini","Carlia_dogare")),
                ln.offset=1.06,lab.offset=1.09, mark.node=FALSE)

#Label padbearing species on tips
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")
#Label padbearing species on tips
Full.Relaxed.Dataset <- read.csv("FullRelaxed_TraitAnalysis.csv",row.names=1)
Full.toepad.tips<-read.csv("FullRelaxed_TraitAnalysis.csv")
Full.toepad.tips <- filter(Full.toepad.tips, toepads == "1")
Full.Toepad.Species <- Full.toepad.tips$Concatenated
ii<-sapply(Full.Toepad.Species,grep,FullTree$tip.label)
ii
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")

#Examine the distribution of estimated origins and losses of arboreality in the Relaxed arboreality dataset

Full.Relaxed.Arboreal.ARD.HPD<-describe.simmap(Full.Relaxed.Arboreal.simmap.ARD)
Full.Relaxed.Arboreal.ER.HPD<-describe.simmap(Full.Relaxed.Arboreal.simmap.ER)

#Summary statistics

Full.Relaxed.Arboreal.simmap.ARD.density <- density(Full.Relaxed.Arboreal.simmap.ARD) #Summary statistics
Full.Relaxed.Arboreal.simmap.ER.density <- density(Full.Relaxed.Arboreal.simmap.ER) #Summary statistics

#Relaxed arboreality frequency distribution
Full.Relaxed.Arbor.ARD.origins.losses<- Full.Relaxed.Arboreal.ARD.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Full.Relaxed.Arbor.ARD.origins.losses) <- c("Gains", "Losses") #rename columns
Full.Relaxed.Arbor.ARD.origins.losses.df <- as.data.frame(Full.Relaxed.Arbor.ARD.origins.losses)

#gains of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)

Gains.Distribution.Full.Relaxed.Arbor.ARD <- ggplot(Full.Relaxed.Arbor.ARD.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_classic()+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+ggtitle("Gains.Distribution.Full.Relaxed.Arbor.ARD")+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 243, xend = 287, 
               y = 23.5, yend = 23.5,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320)) 

#losses of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Full.Relaxed.Arbor.ARD <- ggplot(Full.Relaxed.Arbor.ARD.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_classic()+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+ggtitle("Losses.Distribution.Full.Relaxed.Arbor.ARD")+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 209, xend = 264, 
               y = 23.5, yend = 23.5,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320)) 

#Put both histograms together
Full.Relaxed.Arbor.ARD.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Relaxed.Arbor.ARD, Losses.Distribution.Full.Relaxed.Arbor.ARD,labels = c("A", "B", ncol = 1, nrow = 1))

Full.Relaxed.Arboreal.simmap.ER.density$hpd
#Relaxed arboreality frequency distribution
Full.Relaxed.Arbor.ER.origins.losses<- Full.Relaxed.Arboreal.ER.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Full.Relaxed.Arbor.ER.origins.losses) <- c("Gains", "Losses") #rename columns
Full.Relaxed.Arbor.ER.origins.losses.df <- as.data.frame(Full.Relaxed.Arbor.ER.origins.losses)
#gains of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Full.Relaxed.Arbor.ER <- ggplot(Full.Relaxed.Arbor.ER.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_classic()+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+ggtitle("Gains.Distribution.Full.Relaxed.Arbor.ER")+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 264, xend = 311, 
               y = 24, yend = 24,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320)) 
#losses of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Full.Relaxed.Arbor.ER <- ggplot(Full.Relaxed.Arbor.ER.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_classic()+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+ggtitle("Losses.Distribution.Full.Relaxed.Arbor.ER")+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 174, xend = 224, 
               y = 19.5, yend = 19.5,
               colour = "black")+scale_y_continuous(limits = c(0,20.5), expand = c(0, 0))+coord_cartesian(xlim = c(0, 320))

#Put both histograms together
Full.Relaxed.Arbor.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Relaxed.Arbor.ER, Losses.Distribution.Full.Relaxed.Arbor.ER,labels = c("A", "B", ncol = 1, nrow = 1))

## ARD and ER together
Full.Relaxed.Arbor.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Relaxed.Arbor.ARD, Losses.Distribution.Full.Relaxed.Arbor.ARD, Gains.Distribution.Full.Relaxed.Arbor.ER, Losses.Distribution.Full.Relaxed.Arbor.ER,labels = c("A", "B", "C", "D",ncol = 2, nrow = 2))
Full.Strict.Arbor.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Strict.Arbor.ARD, Losses.Distribution.Full.Strict.Arbor.ARD, Gains.Distribution.Full.Strict.Arbor.ER, Losses.Distribution.Full.Strict.Arbor.ER,labels = c("A", "B", "C", "D", ncol = 2, nrow = 2))
Full.Strict.Relaxed.Arbor.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Strict.Arbor.ARD, Losses.Distribution.Full.Strict.Arbor.ARD, Gains.Distribution.Full.Strict.Arbor.ER, Losses.Distribution.Full.Strict.Arbor.ER,Gains.Distribution.Full.Relaxed.Arbor.ARD, Losses.Distribution.Full.Relaxed.Arbor.ARD, Gains.Distribution.Full.Relaxed.Arbor.ER, Losses.Distribution.Full.Relaxed.Arbor.ER,labels = c("A", "B", "C", "D", "E", "F", "G", "H",ncol = 1, nrow = 2))

grid.arrange(arrangeGrob(Gains.Distribution.Full.Strict.Arbor.ARD, 
                         Losses.Distribution.Full.Strict.Arbor.ARD, 
                         Gains.Distribution.Full.Strict.Arbor.ER, 
                         Losses.Distribution.Full.Strict.Arbor.ER,
                         Gains.Distribution.Full.Relaxed.Arbor.ARD, 
                         Losses.Distribution.Full.Relaxed.Arbor.ARD, 
                         Gains.Distribution.Full.Relaxed.Arbor.ER, 
                         Losses.Distribution.Full.Relaxed.Arbor.ER, ncol = 2), # Second row with 2 plots in 2 different columns
             nrow = 1)  


##########################################################################
######                                                              ######
######   	                  Full Multistate                         ######
######                                                              ######
##########################################################################

Full.Multi.Dataset <- read.csv("Full_Dataset_2692_SupplementalTable.csv",row.names=1)

#Subset some characters for subsequent analyses

Full.Multi.Arboreal <- as.factor(setNames(Full.Multi.Dataset[,6],rownames(Full.Multi.Dataset))) #toepads ful ldataset

#Assign prior probability on root node for non-arboreality

Multi.pi<-setNames(c(0,1,0),c("arboreal", "not.arboreal","semi.arboreal"))

#Multi.pi<-setNames(c(1,0,0),c("not.arboreal", "semi.arboreal","arboreal"))
#fit two Mk models

#fitARD.Multi.arbor.full<-fitMk(FullTree,Full.Multi.Arboreal,pi=Multi.pi,model="ARD")
#AIC(fitARD.Multi.arbor.full)
#
#fitER.Multi.arbor.full<<-fitMk(FullTree,Full.Multi.Arboreal,pi=Multi.pi,model="ER")
#AIC(fitER.Multi.arbor.full)
#
##Compare AIC values and determine best model
#
#aic.Multi.arbor.full<-setNames(sapply(list(fitARD.Multi.arbor.full,fitER.Multi.arbor.full),AIC), c("fitARD","fitER"))
#aic.Multi.arbor.full
#aic.w(aic.Multi.arbor.full)

## Full Multi Arboreal reconstructions

Full.Multi.Arboreal.simmap.ARD<-make.simmap(FullTree,Full.Multi.Arboreal,nsim=500,pi=Multi.pi,model="ARD")
Full.Multi.Arboreal.simmap.ER<-make.simmap(FullTree,Full.Multi.Arboreal,nsim=500,pi=Multi.pi,model="ER")
Full.Multi.Arboreal.simmap.SYM<-make.simmap(FullTree,Full.Multi.Arboreal,nsim=500,pi=Multi.pi,model="SYM")

#plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Full.Multi.Arboreal.simmap.ARD<-summary(Full.Multi.Arboreal.simmap.ARD)
cols<-setNames(c("#D8B365", "#F5F5F5", "#5AB4AC"),levels(Full.Multi.Arboreal))
plot(Obj.Full.Multi.Arboreal.simmap.ARD,fsize=0.6,ftype="off",colors=cols,ylim=c(-2,Ntip(FullTree)))
add.simmap.legend(colors=cols[1:3],prompt=FALSE,x=0,y=2500,vertical=TRUE)
title(main = "Full Multi Arbor ARD 500")
cladelabels(tree=FullTree, "Gekkota", node=2701, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=FullTree, "Dactyloidae", node=5190, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=FullTree, "Scincidae", node=3426, wing.length=NULL, cex=1,offset=1.5)
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey", offset = 1.4)

Obj.Full.Multi.Arboreal.simmap.ER<-summary(Full.Multi.Arboreal.simmap.ER)
cols<-setNames(c("#D8B365", "#F5F5F5", "#5AB4AC"),levels(Full.Multi.Arboreal))
plot(Obj.Full.Multi.Arboreal.simmap.ER,fsize=0.6,ftype="off",colors=cols,ylim=c(-2,Ntip(FullTree)))
add.simmap.legend(colors=cols[1:3],prompt=FALSE,x=0,y=2500,vertical=FALSE)
title(main = "Full Multi Arbor ER 500")
cladelabels(tree=FullTree, "Gekkota", node=2701, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=FullTree, "Dactyloidae", node=5190, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=FullTree, "Scincidae", node=3426, wing.length=NULL, cex=1,offset=1.5)
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey", offset = 1.4)

Obj.Full.Multi.Arboreal.simmap.SYM<-summary(Full.Multi.Arboreal.simmap.SYM)
cols<-setNames(c("#D8B365", "#F5F5F5", "#5AB4AC"),levels(Full.Multi.Arboreal))
plot(Obj.Full.Multi.Arboreal.simmap.SYM,fsize=0.6,ftype="off",colors=cols,ylim=c(-2,Ntip(FullTree)))
add.simmap.legend(colors=cols[1:3],prompt=FALSE,x=0,y=2500,vertical=TRUE)
title(main = "Full Multi Arbor SYM 500")
cladelabels(tree=FullTree, "Gekkota", node=2701, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=FullTree, "Dactyloidae", node=5190, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=FullTree, "Scincidae", node=3426, wing.length=NULL, cex=1,offset=1.5)
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey", offset = 1.4)


#Label padbearing species on tips
Full.toepad.tips<-read.csv("Full_Dataset_2692_SupplementalTable.csv")
Full.toepad.tips <- filter(Full.toepad.tips, toepads == "1")
Full.Toepad.Species <- Full.toepad.tips$Concatenated
ii<-sapply(Full.Toepad.Species,grep,FullTree$tip.label)
ii
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")

#Examine the distribution of estimated origins and losses of arboreality in the Multi arboreality dataset

Full.Multi.Arboreal.ARD.HPD<-describe.simmap(Full.Multi.Arboreal.simmap.ARD)
Full.Multi.Arboreal.ER.HPD<-describe.simmap(Full.Multi.Arboreal.simmap.ER)
Full.Multi.Arboreal.SYM.HPD<-describe.simmap(Full.Multi.Arboreal.simmap.SYM)

Multi.arboreal.HPD
#Summary statistics
Full.Multi.Arboreal.simmap.ARD.density <- density(Full.Multi.Arboreal.simmap.ARD) #Summary statistics
Full.Multi.Arboreal.simmap.ER.density <- density(Full.Multi.Arboreal.simmap.ER) #Summary statistics
Full.Multi.Arboreal.simmap.SYM.density <- density(Full.Multi.Arboreal.simmap.SYM) #Summary statistics



## SYM multi distrib full
Full.Multi.Arbor.SYM.origins.losses <- Full.Multi.Arboreal.SYM.HPD$count[,2:7] #all 500 reconstructions with number of origins and reconstructions for each
colnames(Full.Multi.Arbor.SYM.origins.losses) <- c("Arboreal.to.Non.arboreal", "Arboreal.to.Semi.arboreal","Non.arboreal.to.Arboreal","Non.arboreal.to.Semi.arboreal","Semi.arboreal.to.Arboreal","Semi.arboreal.to.Non.arboreal") #rename columns
Full.Multi.Arbor.SYM.origins.losses <- as.data.frame(Full.Multi.Arbor.SYM.origins.losses)
summary(Full.Multi.Arbor.SYM.origins.losses) #Summary statistics for reconstructions
#HPD Calcs
Arboreal.to.Semi.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Full.Multi.Arbor.SYM.origins.losses))[2,]
#157 - 211
Non.arboreal.to.Semi.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Full.Multi.Arbor.SYM.origins.losses))[4,]
#223 - 280
Semi.arboreal.to.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Full.Multi.Arbor.SYM.origins.losses))[5,]
#154 - 193
Semi.arboreal.to.Non.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Full.Multi.Arbor.SYM.origins.losses))[6,]
#326 - 397
Arboreal.to.Non.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Full.Multi.Arbor.SYM.origins.losses))[1,]
#0-3
#Arboreal to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.SYM.HPD.Arboreal.to.Semi.Arboreal <- ggplot(Full.Multi.Arbor.SYM.origins.losses, aes(x=Arboreal.to.Semi.arboreal))+
  geom_histogram(color="#A6611A", fill="#A6611A",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Semi.arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.SYM")+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 111, xend = 153, 
               y = 24, yend = 24,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 360))+ theme(plot.title = element_text(size=10))
#Non-arboreality to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.SYM.HPD.NonArboreal.to.Semi.Arboreal <- ggplot(Full.Multi.Arbor.SYM.origins.losses, aes(x=Non.arboreal.to.Semi.arboreal))+
  geom_histogram(color="#DFC27D", size=0.2,fill="#DFC27D",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Semi.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.SYM")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 294, xend = 346, 
               y = 23, yend = 23,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 360))+ theme(plot.title = element_text(size=10))


#Semi-arboreality to arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.SYM.HPD.SemiArboreal.to.Arboreal <- ggplot(Full.Multi.Arbor.SYM.origins.losses, aes(x=Semi.arboreal.to.Arboreal))+
  geom_histogram(color="gray", size=0.2,fill="gray",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.SYM")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 185, xend = 225, 
               y = 28, yend = 28,
               colour = "black")+scale_y_continuous(limits = c(0,32), expand = c(0, 0))+coord_cartesian(xlim = c(0, 360))+ theme(plot.title = element_text(size=10))

#Semi-arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.SYM.HPD.SemiArboreal.to.Non.Arboreal <- ggplot(Full.Multi.Arbor.SYM.origins.losses, aes(x=Semi.arboreal.to.Non.arboreal))+
  geom_histogram(color="#80CDC1", size=0.2,fill="#80CDC1",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.SYM")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 153, xend = 202, 
               y = 23.5, yend = 23.5,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 360))+ theme(plot.title = element_text(size=10))


#Put everything together
Multistate.SYM.Full <- grid.arrange(Full.Multi.SYM.HPD.SemiArboreal.to.Non.Arboreal,
                                    Full.Multi.SYM.HPD.SemiArboreal.to.Arboreal,
                                    Full.Multi.SYM.HPD.NonArboreal.to.Semi.Arboreal,
                                    Full.Multi.SYM.HPD.Arboreal.to.Semi.Arboreal, layout_matrix = rbind(c(1,1,1),c(2,3,4)))

Multistate.SYM.Full <- ggarrange(Full.Multi.SYM.HPD.SemiArboreal.to.Non.Arboreal,
                                Full.Multi.SYM.HPD.SemiArboreal.to.Arboreal,
                                Full.Multi.SYM.HPD.NonArboreal.to.Semi.Arboreal,
                                Full.Multi.SYM.HPD.Arboreal.to.Semi.Arboreal,labels = c("A", "B", "C", "D", ncol = 3, nrow = 3))

## ER multi distrib full
Full.Multi.Arbor.ER.origins.losses <- Full.Multi.Arboreal.ER.HPD$count[,2:7] #all 500 reconstructions with number of origins and reconstructions for each
colnames(Full.Multi.Arbor.ER.origins.losses) <- c("Arboreal.to.Non.arboreal", "Arboreal.to.Semi.arboreal","Non.arboreal.to.Arboreal","Non.arboreal.to.Semi.arboreal","Semi.arboreal.to.Arboreal","Semi.arboreal.to.Non.arboreal") #rename columns
Full.Multi.Arbor.ER.origins.losses <- as.data.frame(Full.Multi.Arbor.ER.origins.losses)
summary(Full.Multi.Arbor.ER.origins.losses) #Summary statistics for reconstructions
#HPD Calcs
Arboreal.to.Semi.arboreal.HPD.ER <- HPDinterval(as.mcmc(Full.Multi.Arbor.ER.origins.losses))[2,]
#157 - 211
Non.arboreal.to.Semi.arboreal.HPD.ER <- HPDinterval(as.mcmc(Full.Multi.Arbor.ER.origins.losses))[4,]
#223 - 280
Semi.arboreal.to.arboreal.HPD.ER <- HPDinterval(as.mcmc(Full.Multi.Arbor.ER.origins.losses))[5,]
#154 - 193
Semi.arboreal.to.Non.arboreal.HPD.ER <- HPDinterval(as.mcmc(Full.Multi.Arbor.ER.origins.losses))[6,]
#326 - 397
Arboreal.to.Non.arboreal.HPD.ER <- HPDinterval(as.mcmc(Full.Multi.Arbor.ER.origins.losses))[1,]
#0-3
#Arboreal to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ER.HPD.Arboreal.to.Semi.Arboreal <- ggplot(Full.Multi.Arbor.ER.origins.losses, aes(x=Arboreal.to.Semi.arboreal))+
  geom_histogram(color="#8C510A", fill="#8C510A",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Semi.arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ER")+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 80, xend = 107, 
               y = 32, yend = 32,
               colour = "black")+scale_y_continuous(limits = c(0,35), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))
#Non-arboreality to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ER.HPD.NonArboreal.to.Semi.Arboreal <- ggplot(Full.Multi.Arbor.ER.origins.losses, aes(x=Non.arboreal.to.Semi.arboreal))+
  geom_histogram(color="#D8B365", size=0.2,fill="#D8B365",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Semi.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 241, xend = 281, 
               y = 26, yend = 26,
               colour = "black")+scale_y_continuous(limits = c(0,27), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))


#Semi-arboreality to arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ER.HPD.SemiArboreal.to.Arboreal <- ggplot(Full.Multi.Arbor.ER.origins.losses, aes(x=Semi.arboreal.to.Arboreal))+
  geom_histogram(color="gray", size=0.2,fill="gray",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 64, xend = 92, 
               y = 36, yend = 36,
               colour = "black")+scale_y_continuous(limits = c(0,38), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))

#Semi-arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ER.HPD.SemiArboreal.to.Non.Arboreal <- ggplot(Full.Multi.Arbor.ER.origins.losses, aes(x=Semi.arboreal.to.Non.arboreal))+
  geom_histogram(color="#C7EAE5", size=0.2,fill="#C7EAE5",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 67, xend = 104, 
               y = 33, yend = 33,
               colour = "black")+scale_y_continuous(limits = c(0,35), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))

#Arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ER.HPD.Arboreal.to.Non.Arboreal <- ggplot(Full.Multi.Arbor.ER.origins.losses, aes(x=Arboreal.to.Non.arboreal))+
  geom_histogram(color="#5AB4AC", size=0.2,fill="#5AB4AC",binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Non.arboreal)),
                                                                                                             color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 43, xend = 68, 
               y = 37, yend = 37,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))

#Non-arboreality to Arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ER.HPD.Non.Arboreal.to.Arboreal <- ggplot(Full.Multi.Arbor.ER.origins.losses, aes(x=Non.arboreal.to.Arboreal))+
  geom_histogram(color="#01665E", size=0.2,fill="#01665E",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Arboreal)),
                                                                                                             color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 82, xend = 110, 
               y = 34, yend = 34,
               colour = "black")+scale_y_continuous(limits = c(0,35), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))


#Put everything together
#Multistate.ER.Full <- grid.arrange(Full.Multi.ER.HPD.SemiArboreal.to.Non.Arboreal,
#                                    Full.Multi.ER.HPD.Arboreal.to.Non.Arboreal, 
#                                    Full.Multi.ER.HPD.SemiArboreal.to.Arboreal,
#                                    Full.Multi.ER.HPD.NonArboreal.to.Semi.Arboreal,
#                                    Full.Multi.ER.HPD.Arboreal.to.Semi.Arboreal, 
#                                    Full.Multi.ER.HPD.Non.Arboreal.to.Arboreal,nrow=3,ncol=3)
Multistate.ER.Full <- ggarrange(Full.Multi.ER.HPD.SemiArboreal.to.Non.Arboreal,
                                                    Full.Multi.ER.HPD.Arboreal.to.Non.Arboreal, 
                                                    Full.Multi.ER.HPD.SemiArboreal.to.Arboreal,
                                                    Full.Multi.ER.HPD.NonArboreal.to.Semi.Arboreal,
                                                    Full.Multi.ER.HPD.Arboreal.to.Semi.Arboreal, 
                                                    Full.Multi.ER.HPD.Non.Arboreal.to.Arboreal,labels = c("A", "B", "C", "D", "E", "F", ncol = 3, nrow = 3))



## ARD multi distrib full
Full.Multi.Arbor.ARD.origins.losses <- Full.Multi.Arboreal.ARD.HPD$count[,2:7] #all 500 reconstructions with number of origins and reconstructions for each
colnames(Full.Multi.Arbor.ARD.origins.losses) <- c("Arboreal.to.Non.arboreal", "Arboreal.to.Semi.arboreal","Non.arboreal.to.Arboreal","Non.arboreal.to.Semi.arboreal","Semi.arboreal.to.Arboreal","Semi.arboreal.to.Non.arboreal") #rename columns
Full.Multi.Arbor.ARD.origins.losses <- as.data.frame(Full.Multi.Arbor.ARD.origins.losses)
summary(Full.Multi.Arbor.ARD.origins.losses) #Summary statistics for reconstructions
#HPD Calcs
Arboreal.to.Semi.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Full.Multi.Arbor.ARD.origins.losses))[2,]
#157 - 211
Non.arboreal.to.Semi.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Full.Multi.Arbor.ARD.origins.losses))[4,]
#223 - 280
Semi.arboreal.to.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Full.Multi.Arbor.ARD.origins.losses))[5,]
#154 - 193
Semi.arboreal.to.Non.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Full.Multi.Arbor.ARD.origins.losses))[6,]
#326 - 397
Arboreal.to.Non.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Full.Multi.Arbor.ARD.origins.losses))[1,]
#0-3
#Arboreal to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ARD.HPD.Arboreal.to.Semi.Arboreal <- ggplot(Full.Multi.Arbor.ARD.origins.losses, aes(x=Arboreal.to.Semi.arboreal))+
  geom_histogram(color="#A6611A", fill="#A6611A",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Semi.arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ARD")+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 157, xend = 211, 
               y = 23, yend = 23,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))
#Non-arboreality to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ARD.HPD.NonArboreal.to.Semi.Arboreal <- ggplot(Full.Multi.Arbor.ARD.origins.losses, aes(x=Non.arboreal.to.Semi.arboreal))+
  geom_histogram(color="#DFC27D", size=0.2,fill="#DFC27D",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Semi.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 223, xend = 280, 
               y = 24.5, yend = 24.5,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))


#Semi-arboreality to arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ARD.HPD.SemiArboreal.to.Arboreal <- ggplot(Full.Multi.Arbor.ARD.origins.losses, aes(x=Semi.arboreal.to.Arboreal))+
  geom_histogram(color="gray", size=0.2,fill="gray",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 154, xend = 193, 
               y = 35, yend = 35,
               colour = "black")+scale_y_continuous(limits = c(0,36), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))

#Semi-arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ARD.HPD.SemiArboreal.to.Non.Arboreal <- ggplot(Full.Multi.Arbor.ARD.origins.losses, aes(x=Semi.arboreal.to.Non.arboreal))+
  geom_histogram(color="#80CDC1", size=0.2,fill="#80CDC1",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 326, xend = 397, 
               y = 17, yend = 17,
               colour = "black")+scale_y_continuous(limits = c(0,18), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))

#Arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Full.Multi.ARD.HPD.Arboreal.to.Non.Arboreal <- ggplot(Full.Multi.Arbor.ARD.origins.losses, aes(x=Arboreal.to.Non.arboreal))+
  geom_histogram(color="#018571", size=0.2,fill="#018571",binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Full.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 0, xend = 3, 
               y = 197, yend = 197,
               colour = "black")+scale_y_continuous(limits = c(0,200), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))



#Put everything together
Multistate.ARD.Full <- grid.arrange(Full.Multi.ARD.HPD.SemiArboreal.to.Non.Arboreal,
             Full.Multi.ARD.HPD.Arboreal.to.Non.Arboreal, 
             Full.Multi.ARD.HPD.SemiArboreal.to.Arboreal,
             Full.Multi.ARD.HPD.NonArboreal.to.Semi.Arboreal,
             Full.Multi.ARD.HPD.Arboreal.to.Semi.Arboreal, layout_matrix = rbind(c(1,1,1,1),c(2,3,4,5)))



##########################################################################
######                                                              ######
######   	                  Full Toepads                            ######
######                                                              ######
##########################################################################

Full.Toepads <- as.factor(setNames(Full.Strict.Dataset[,2],rownames(Full.Strict.Dataset))) #toepads ful ldataset

#Assign prior probability on root node for non-arboreality

toepads.pi<-setNames(c(1,0),c("0","1"))

#fit two Mk models

#fitARD.Toepads.full<-fitMk(FullTree,Full.Toepads,pi=toepads.pi,model="ARD")
#AIC(fitARD.Toepads.full)
#
#fitER.Toepads.full<<-fitMk(FullTree,Full.Toepads,pi=toepads.pi,model="ER")
#AIC(fitER.Toepads.full)
#
##Compare AIC values and determine best model
#
#aic.Toepads.full<-setNames(sapply(list(fitARD.Toepads.full,fitER.Toepads.full),AIC), c("fitARD","fitER"))
#aic.Toepads.full
#aic.w(aic.Toepads.full)

## Full Strict Arboreal reconstructions

Full.Toepads.simmap.ARD<-make.simmap(FullTree,Full.Toepads,nsim=500,pi=toepads.pi,model="ARD")
Full.Toepads.simmap.ER<-make.simmap(FullTree,Full.Toepads,nsim=500,pi=toepads.pi,model="ER")

#plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Full.Toepads.simmap.ARD <-densityMap(Full.Toepads.simmap.ARD,states=levels(Full.Toepads)[1:2],plot=FALSE)
colors.5<-setMap(Obj.Full.Toepads.simmap.ARD,c("gray84","lightskyblue"))
plot(colors.5,Obj.Full.Toepads.simmap.ARD,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Full Toepads ARD 500")

###Major padbearing clades
#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(FullTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=findMRCA(FullTree,c("Typhlosaurus_braini","Carlia_dogare")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(FullTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)


Obj.Full.Toepads.simmap.ER <-densityMap(Full.Toepads.simmap.ER,states=levels(Full.Toepads)[1:2],plot=FALSE)
colors.6<-setMap(Obj.Full.Toepads.simmap.ER,c("gray84","lightskyblue"))
plot(colors.6,Obj.Full.Toepads.simmap.ER,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Full Toepads ER 500")

#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(FullTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(FullTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=findMRCA(FullTree,c("Typhlosaurus_braini","Carlia_dogare")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Label strict arboreal species on tips
#Full.Strict.Arboreal.tips<-read.csv("Full_Dataset_2692_SupplementalTable.csv")
#Full.Strict.Arboreal.tips <- filter(Full.Strict.Arboreal.tips, strict.arboreal == "1")
#Full.Strict.Arboreal.tips.Species <- Full.Strict.Arboreal.tips$Concatenated
#arbor.ii<-sapply(Full.Strict.Arboreal.tips.Species,grep,FullTree$tip.label)
#arbor.ii
#ape::tiplabels(pch=21, tip=arbor.ii, cex=2, col='black',bg="grey")

#Examine the distribution of estimated origins and losses of arboreality in the strict arboreality dataset

Full.Toepads.ARD.HPD<-describe.simmap(Full.Toepads.simmap.ARD)
Full.Toepads.ER.HPD<-describe.simmap(Full.Toepads.simmap.ER)

#Summary statistics
Full.Toepads.simmap.ARD.density <- density(Full.Toepads.simmap.ARD) #Summary statistics
Full.Toepads.simmap.ER.density <- density(Full.Toepads.simmap.ER) #Summary statistics

#Strict arboreality frequency distribution
Full.Toepads.ARD.origins.losses<- Full.Toepads.ARD.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Full.Toepads.ARD.origins.losses) <- c("Gains", "Losses") #rename columns
Full.Toepads.ARD.origins.losses.df <- as.data.frame(Full.Toepads.ARD.origins.losses)
#gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Full.Toepads.ARD <- ggplot(Full.Toepads.ARD.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Toepads gains") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+ggtitle("Gains.Distribution.Full.Toepads.ARD")+
  geom_text(x = 12,  y = 340, 
            label = "", 
            colour = "black")+
  geom_segment(x = 3, xend = 5, 
               y = 330, yend = 330,
               colour = "black")+scale_y_continuous(limits = c(0,350), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30))
#losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Full.Toepads.ARD <- ggplot(Full.Toepads.ARD.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="black", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+
  xlab("Toepad losses") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+ggtitle("Losses.Distribution.Full.Toepads.ARD")+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 19, xend = 26, 
               y = 140, yend = 140,
               colour = "black")+scale_y_continuous(limits = c(0,150), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30)) 

#Put both histograms together
Full.Toepads.ARD.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Toepads.ARD, Losses.Distribution.Full.Toepads.ARD,labels = c("A", "B", ncol = 1, nrow = 1))

#Strict arboreality frequency distribution
Full.Toepads.ER.origins.losses<- Full.Toepads.ER.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Full.Toepads.ER.origins.losses) <- c("Gains", "Losses") #rename columns
Full.Toepads.ER.origins.losses.df <- as.data.frame(Full.Toepads.ER.origins.losses)
Full.Toepads.simmap.ER.density
#gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Full.Toepads.ER <- ggplot(Full.Toepads.ER.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Toepad gains") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+ggtitle("Gains.Distribution.Full.Toepads.ER")+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 9, xend = 12, 
               y = 320, yend = 320,
               colour = "black")+scale_y_continuous(limits = c(0,340), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30)) 
#losses of toepads with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Full.Toepads.ER <- ggplot(Full.Toepads.ER.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="black", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+
  xlab("Toepad losses") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+ggtitle("Losses.Distribution.Full.Toepads.ER")+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 10, xend = 13, 
               y = 290, yend = 290,
               colour = "black")+scale_y_continuous(limits = c(0,300), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30)) 

#Put both histograms together
Full.Toepads.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Toepads.ER, Losses.Distribution.Full.Toepads.ER,labels = c("A", "B", ncol = 1, nrow = 1))

## Full toepads ARD/ER
Full.Toepads.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Full.Toepads.ARD, Losses.Distribution.Full.Toepads.ARD, Gains.Distribution.Full.Toepads.ER, Losses.Distribution.Full.Toepads.ER, labels = c("A", "B", "C", "D", ncol = 2, nrow = 2))






