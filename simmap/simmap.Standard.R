
##########################################################################
######                                                              ######
######   	            stochastic character mapping                  ######
######                                                              ######
##########################################################################

library(phytools)
library(ape)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(coda)
#setwd

## Standard Dataset

## Read in trait datasets

Standard.Strict.Dataset<-read.csv("StandardStrict_TraitAnalysis.csv",row.names=1)
Standard.Relaxed.Dataset<-read.csv("StandardRelaxed_TraitAnalysis.csv",row.names=1)

## Read in tree 
StandardTree <- read.nexus(file = "StandardTree.tre")

## Subset some characters for subsequent analyses

Standard.Strict.Arboreal <- as.factor(setNames(Standard.Strict.Dataset[,1],rownames(Standard.Strict.Dataset))) #Strict arboreal Standard dataset
Standard.Toepads <- as.factor(setNames(Standard.Strict.Dataset[,2],rownames(Standard.Strict.Dataset))) #toepads Standard dataset

## Assign prior probability on root node for non-arboreality

arbor.pi<-setNames(c(1,0),c("0","1"))
arbor.pi

# fit two Mk models

fitARD.strict.arbor.Standard<-fitMk(StandardTree,Standard.Strict.Arboreal,pi=arbor.pi,model="ARD")
AIC(fitARD.strict.arbor.Standard)

fitER.strict.arbor.Standard<<-fitMk(StandardTree,Standard.Strict.Arboreal,pi=arbor.pi,model="ER")
AIC(fitER.strict.arbor.Standard)

# Compare AIC values and determine best model

aic.strict.arbor.Standard<-setNames(sapply(list(fitARD.strict.arbor.Standard,fitER.strict.arbor.Standard),AIC), c("fitARD","fitER"))
aic.strict.arbor.Standard
aic.w(aic.strict.arbor.Standard)

# Standard Strict Arboreal reconstructions
Standard.Strict.Arboreal.simmap.ARD<-make.simmap(StandardTree,Standard.Strict.Arboreal,nsim=500,pi=arbor.pi,model="ARD")
Standard.Strict.Arboreal.simmap.ER<-make.simmap(StandardTree,Standard.Strict.Arboreal,nsim=500,pi=arbor.pi,model="ER")

## plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Standard.Strict.Arboreal.simmap.ARD <-densityMap(Standard.Strict.Arboreal.simmap.ARD,states=levels(Standard.Strict.Arboreal)[1:2],plot=FALSE)
colors.1<-setMap(Obj.Standard.Strict.Arboreal.simmap.ARD,c("gray84","lightskyblue"))
plot(colors.1,Obj.Standard.Strict.Arboreal.simmap.ARD,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Standard.Strict.Arboreal.simmap.ARD")
ape::tiplabels(pch=21, tip=ii, cex=1.5, col='black',bg="grey")

#Plot some family arc labels

#Agamidae
#arc.cladelabels(text="Agamidae",node=findMRCA(StandardTree,c("Saara_asmussi","Pseudocalotes_brevipes")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Anguidae
#arc.cladelabels(text="Anguidae",node=findMRCA(StandardTree,c("Abronia_campbelli","Celestus_enneagrammus")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Chamaeleonidae
#arc.cladelabels(text="Chamaeleonidae",node=findMRCA(StandardTree,c("Trioceros_rudis","Palleon_nasus")),
#                ln.offset=1.05,lab.offset=1.07,orientation = "horizontal", mark.node=FALSE)
##Corytophanidae
#arc.cladelabels(text="Corytophanidae",node=findMRCA(StandardTree,c("Basiliscus_plumifrons","Laemanctus_longipes")),
#                ln.offset=1.09,lab.offset=1.12, orientation = "horizontal",mark.node=FALSE)
##Dactyloidae
#arc.cladelabels(text="Dactyloidae",node=findMRCA(StandardTree,c("Anolis_bonairensis","Anolis_poecilopus")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Diplodactylidae
#arc.cladelabels(text="Diplodactylidae",node=findMRCA(StandardTree,c("Pseudothecadactylus_lindneri","Diplodactylus_pulcher")),
#                ln.offset=1.03,lab.offset=1.06,mark.node=FALSE)
##Eublepharidae
#arc.cladelabels(text="Eublepharidae",node=findMRCA(StandardTree,c("Aeluroscalabotes_felinus","Hemitheconyx_caudicinctus")),
#                ln.offset=1.15,lab.offset=1.18,orientation = "horizontal", mark.node=FALSE)
##Gekkonidae
#arc.cladelabels(text="Gekkonidae",node=findMRCA(StandardTree,c("Lepidodactylus_novaeguineae","Phelsuma_lineata")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Gerrhosauridae
#arc.cladelabels(text="Gerrhosauridae",node=findMRCA(StandardTree,c("Broadleysaurus_major","Zonosaurus_bemaraha")),
#                ln.offset=1.03,lab.offset=1.05,orientation = "horizontal", mark.node=FALSE)
##Iguanidae
#arc.cladelabels(text="Iguanidae",node=findMRCA(StandardTree,c("Dipsosaurus_dorsalis","Ctenosaura_palearis")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Lacertidae
#arc.cladelabels(text="Lacertidae",node=findMRCA(StandardTree,c("Psammodromus_hispanicus","Darevskia_caucasica")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Leiosauridae
#arc.cladelabels(text="Leiosauridae",node=findMRCA(StandardTree,c("Enyalius_bilineatus","Leiosaurus_paronae")),
#                ln.offset=1.21,lab.offset=1.24, orientation = "horizontal",mark.node=FALSE)
#
##Liolaemidae
#arc.cladelabels(text="Liolaemidae",node=findMRCA(StandardTree,c("Ctenoblepharys_adspersa","Liolaemus_ornatus")),
#                ln.offset=1.03,lab.offset=1.06, orientation = "horizontal",mark.node=FALSE)
#
### PICK UP HERE
##Opluridae
#arc.cladelabels(text="Opluridae",node=findMRCA(StandardTree,c("Chalarodon_madagascariensis","Oplurus_grandidieri")),
#                ln.offset=1.15,lab.offset=1.18,orientation = "horizontal", mark.node=FALSE)
##Phrynosomatidae
#arc.cladelabels(text="Phrynosomatidae",node=findMRCA(StandardTree,c("Callisaurus_draconoides","Sceloporus_stejnegeri")),
#                ln.offset=1.09,lab.offset=1.12,orientation = "horizontal", mark.node=FALSE)
##Phyllodactylidae
#arc.cladelabels(text="Phyllodactylidae",node=findMRCA(StandardTree,c("Thecadactylus_rapicauda","Tarentola_gigas")),
#                ln.offset=1.05,lab.offset=1.08,orientation = "horizontal", mark.node=FALSE)
##Polychrotidae
#arc.cladelabels(text="Polychrotidae",node=findMRCA(StandardTree,c("Polychrus_gutturosus","Polychrus_marmoratus")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Scincidae
#arc.cladelabels(text="Scincidae",node=3052, ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Sphaerodactylidae
#arc.cladelabels(text="Sphaerodactylidae",node=findMRCA(StandardTree,c("Sphaerodactylus_ocoae","Quedenfeldtia_moerens")),
#                ln.offset=1.03,lab.offset=1.06, orientation = "horizontal", mark.node=FALSE)
##Teiidae
#arc.cladelabels(text="Teiidae",node=findMRCA(StandardTree,c("Pholidoscelis_fuscatus","Callopistes_maculatus")),
#                ln.offset=1.05,lab.offset=1.08, mark.node=FALSE)
##Tropiduridae
#arc.cladelabels(text="Tropiduridae",node=findMRCA(StandardTree,c("Stenocercus_ochoai","Tropidurus_cocorobensis")),
#                ln.offset=1.05,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Varanidae
#arc.cladelabels(text="Varanidae",node=findMRCA(StandardTree,c("Varanus_griseus","Varanus_acanthurus")),
#                ln.offset=1.09,lab.offset=1.12, mark.node=FALSE)
#

#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(StandardTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(StandardTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=3052, ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)


Obj.Standard.Strict.Arboreal.simmap.ER <-densityMap(Standard.Strict.Arboreal.simmap.ER,states=levels(Standard.Strict.Arboreal)[1:2],plot=FALSE)
colors.2<-setMap(Obj.Standard.Strict.Arboreal.simmap.ER,c("gray84","lightskyblue"))
plot(colors.2,Obj.Standard.Strict.Arboreal.simmap.ER,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Standard.Strict.Arboreal.simmap.ER")

## Label padbearing species on tips
Standard.toepad.tips<-read.csv("StandardStrict_TraitAnalysis.csv")
Standard.toepad.tips <- filter(Standard.toepad.tips, toepads == "1")
Standard.Toepad.Species <- Standard.toepad.tips$Concatenated
ii<-sapply(Standard.Toepad.Species,grep,StandardTree$tip.label)
ii
ape::tiplabels(pch=21, tip=ii, cex=1, col='black',bg="grey")

#Agamidae
#arc.cladelabels(text="Agamidae",node=findMRCA(StandardTree,c("Uromastyx_princeps","Pseudocalotes_brevipes")),
#                ln.offset=1.02,lab.offset=1.05,orientation = "horizontal", mark.node=FALSE)
##Anguidae
#arc.cladelabels(text="Anguidae",node=findMRCA(StandardTree,c("Abronia_campbelli","Anniella_pulchra")),
#                ln.offset=1.06,lab.offset=1.09,orientation = "horizontal", mark.node=FALSE)
##Chamaeleonidae
#arc.cladelabels(text="Chamaeleonidae",node=findMRCA(StandardTree,c("Chamaeleo_calcaricarens","Brookesia_dentata")),
#                ln.offset=1.07,lab.offset=1.10,orientation = "horizontal", mark.node=FALSE)
##Corytophanidae
#arc.cladelabels(text="Corytophanidae",node=findMRCA(StandardTree,c("Basiliscus_plumifrons","Laemanctus_longipes")),
#                ln.offset=1.025,lab.offset=1.10, orientation = "horizontal",mark.node=FALSE)
##Dactyloidae
#arc.cladelabels(text="Dactyloidae",node=findMRCA(StandardTree,c("Anolis_bonairensis","Anolis_trachyderma")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Diplodactylidae
#arc.cladelabels(text="Diplodactylidae",node=findMRCA(StandardTree,c("Pseudothecadactylus_lindneri","Diplodactylus_pulcher")),
#                ln.offset=1.02,lab.offset=1.05,orientation = "horizontal", mark.node=FALSE)
##Eublepharidae
#arc.cladelabels(text="Eublepharidae",node=findMRCA(StandardTree,c("Aeluroscalabotes_felinus","Hemitheconyx_caudicinctus")),
#                ln.offset=1.08,lab.offset=1.11,orientation = "horizontal", mark.node=FALSE)
##Gekkonidae
#arc.cladelabels(text="Gekkonidae",node=findMRCA(StandardTree,c("Lepidodactylus_novaeguineae","Phelsuma_lineata")),
#                ln.offset=1.03,lab.offset=1.06,orientation = "horizontal", mark.node=FALSE)
##Gerrhosauridae
#arc.cladelabels(text="Gerrhosauridae",node=findMRCA(StandardTree,c("Cordylosaurus_subtessellatus","Zonosaurus_bemaraha")),
#                ln.offset=1.06,lab.offset=1.09,orientation = "horizontal", mark.node=FALSE)
##Iguanidae
#arc.cladelabels(text="Iguanidae",node=findMRCA(StandardTree,c("Dipsosaurus_dorsalis","Ctenosaura_palearis")),
#                ln.offset=1.02,lab.offset=1.05,orientation = "horizontal", mark.node=FALSE)
##Lacertidae
#arc.cladelabels(text="Lacertidae",node=findMRCA(StandardTree,c("Psammodromus_hispanicus","Darevskia_caucasica")),
#                ln.offset=1.02,lab.offset=1.05,orientation = "horizontal", mark.node=FALSE)
##Leiosauridae
#arc.cladelabels(text="Leiosauridae",node=findMRCA(StandardTree,c("Enyalius_bilineatus","Leiosaurus_paronae")),
#                ln.offset=1.06,lab.offset=1.09,orientation = "horizontal", mark.node=FALSE)
##Liolaemidae
#arc.cladelabels(text="Liolaemidae",node=findMRCA(StandardTree,c("Ctenoblepharys_adspersa","Liolaemus_ornatus")),
#                ln.offset=1.06,lab.offset=1.09,orientation = "horizontal", mark.node=FALSE)
##Opluridae
#arc.cladelabels(text="Opluridae",node=findMRCA(StandardTree,c("Chalarodon_madagascariensis","Oplurus_grandidieri")),
#                ln.offset=1.06,lab.offset=1.09, orientation = "horizontal", mark.node=FALSE)
##Phrynosomatidae
#arc.cladelabels(text="Phrynosomatidae",node=findMRCA(StandardTree,c("Callisaurus_draconoides","Sceloporus_stejnegeri")),
#                ln.offset=1.10,lab.offset=1.13, mark.node=FALSE)
##Phyllodactylidae
#arc.cladelabels(text="Phyllodactylidae",node=findMRCA(StandardTree,c("Thecadactylus_rapicauda","Tarentola_gigas")),
#                ln.offset=1.10,lab.offset=1.13, mark.node=FALSE)
##Polychrotidae
#arc.cladelabels(text="Polychrotidae",node=findMRCA(StandardTree,c("Polychrus_gutturosus","Polychrus_marmoratus")),
#                ln.offset=1.10,lab.offset=1.13,orientation = "horizontal", mark.node=FALSE)
##Scincidae
#arc.cladelabels(text="Scincidae",node=3052, ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#
##Sphaerodactylidae
#arc.cladelabels(text="Sphaerodactylidae",node=findMRCA(StandardTree,c("Sphaerodactylus_ocoae","Quedenfeldtia_moerens")),
#                ln.offset=1.02,lab.offset=1.05, mark.node=FALSE)
##Teiidae
#arc.cladelabels(text="Teiidae",node=findMRCA(StandardTree,c("Dicrodon_guttulatum","Callopistes_maculatus")),
#                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
##Tropiduridae
#arc.cladelabels(text="Tropiduridae",node=findMRCA(StandardTree,c("Stenocercus_ochoai","Tropidurus_cocorobensis")),
#                ln.offset=1.02,lab.offset=1.05, mark.node=FALSE)
##Varanidae
#arc.cladelabels(text="Varanidae",node=findMRCA(StandardTree,c("Varanus_griseus","Varanus_acanthurus")),
#                ln.offset=1.02,lab.offset=1.05, mark.node=FALSE)


## Examine the distribution of estimated origins and losses of arboreality in the strict arboreality dataset

Standard.Strict.Arboreal.ARD.HPD<-describe.simmap(Standard.Strict.Arboreal.simmap.ARD)
Standard.Strict.Arboreal.ER.HPD<-describe.simmap(Standard.Strict.Arboreal.simmap.ER)

## Summary statistics
Standard.Strict.Arboreal.simmap.ARD.density <- density(Standard.Strict.Arboreal.simmap.ARD) #Summary statistics
Standard.Strict.Arboreal.simmap.ER.density <- density(Standard.Strict.Arboreal.simmap.ER) #Summary statistics

## Strict arboreality frequency distribution
Standard.Strict.Arbor.ARD.origins.losses<- Standard.Strict.Arboreal.ARD.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Standard.Strict.Arbor.ARD.origins.losses) <- c("Gains", "Losses") #rename columns
Standard.Strict.Arbor.ARD.origins.losses.df <- as.data.frame(Standard.Strict.Arbor.ARD.origins.losses)
# gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Standard.Strict.Arbor.ARD <- ggplot(Standard.Strict.Arbor.ARD.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+ggtitle("Gains.Distribution.Standard.Strict.Arbor.ARD")+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
  label = "", 
  colour = "black")+
  geom_segment(x = 109, xend = 144, 
               y = 33, yend = 33,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))
## losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Standard.Strict.Arbor.ARD <- ggplot(Standard.Strict.Arbor.ARD.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+ggtitle("Losses.Distribution.Standard.Strict.Arbor.ARD")+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 171, xend = 243, 
               y = 19.5, yend = 19.5,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))

## Put both histograms together
Standard.Strict.Arbor.ARD.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Strict.Arbor.ARD, Losses.Distribution.Standard.Strict.Arbor.ARD,labels = c("A", "B", ncol = 1, nrow = 1))


## Strict arboreality frequency distribution
Standard.Strict.Arbor.ER.origins.losses<- Standard.Strict.Arboreal.ER.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Standard.Strict.Arbor.ER.origins.losses) <- c("Gains", "Losses") #rename columns
Standard.Strict.Arbor.ER.origins.losses.df <- as.data.frame(Standard.Strict.Arbor.ER.origins.losses)
## gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Standard.Strict.Arbor.ER <- ggplot(Standard.Strict.Arbor.ER.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+ggtitle("Gains.Distribution.Standard.Strict.Arbor.ER")+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 156, xend = 181, 
               y = 37, yend = 37,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))
## losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Standard.Strict.Arbor.ER <- ggplot(Standard.Strict.Arbor.ER.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+ggtitle("Losses.Distribution.Standard.Strict.Arbor.ER")+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 70, xend = 98, 
               y = 35, yend = 35,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))

## Put both histograms together
Standard.Strict.Arbor.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Strict.Arbor.ER, Losses.Distribution.Standard.Strict.Arbor.ER,labels = c("A", "B", ncol = 1, nrow = 1))
Standard.Strict.Arbor.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Strict.Arbor.ARD, Losses.Distribution.Standard.Strict.Arbor.ARD,Gains.Distribution.Standard.Strict.Arbor.ER, Losses.Distribution.Standard.Strict.Arbor.ER,labels = c("A", "B", "C", "D", ncol = 2, nrow = 2))


##########################################################################
######                                                              ######
######   	          Standard Relaxed Arboreality                    ######
######                                                              ######
##########################################################################

Standard.Relaxed.Dataset <- read.csv("StandardRelaxed_TraitAnalysis.csv",row.names=1)

## Subset some characters for subsequent analyses

Standard.Relaxed.Arboreal <- as.factor(setNames(Standard.Relaxed.Dataset[,1],rownames(Standard.Relaxed.Dataset))) #Relaxed arboreal Standard dataset
Standard.Toepads <- as.factor(setNames(Standard.Relaxed.Dataset[,2],rownames(Standard.Relaxed.Dataset))) #toepads Standard dataset

## Assign prior probability on root node for non-arboreality

arbor.pi<-setNames(c(1,0),c("0","1"))
arbor.pi

## fit two Mk models

#fitARD.Relaxed.arbor.Standard<-fitMk(StandardTree,Standard.Relaxed.Arboreal,pi=arbor.pi,model="ARD")
#AIC(fitARD.Relaxed.arbor.Standard)
#
#fitER.Relaxed.arbor.Standard<<-fitMk(StandardTree,Standard.Relaxed.Arboreal,pi=arbor.pi,model="ER")
#AIC(fitER.Relaxed.arbor.Standard)
#
### Compare AIC values and determine best model
#
#aic.Relaxed.arbor.Standard<-setNames(sapply(list(fitARD.Relaxed.arbor.Standard,fitER.Relaxed.arbor.Standard),AIC), c("fitARD","fitER"))
#aic.Relaxed.arbor.Standard
#aic.w(aic.Relaxed.arbor.Standard)

## Standard Relaxed Arboreal reconstructions

Standard.Relaxed.Arboreal.simmap.ARD<-make.simmap(StandardTree,Standard.Relaxed.Arboreal,nsim=500,pi=arbor.pi,model="ARD")
Standard.Relaxed.Arboreal.simmap.ER<-make.simmap(StandardTree,Standard.Relaxed.Arboreal,nsim=500,pi=arbor.pi,model="ER")

## plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Standard.Relaxed.Arboreal.simmap.ARD <-densityMap(Standard.Relaxed.Arboreal.simmap.ARD,states=levels(Standard.Relaxed.Arboreal)[1:2],plot=FALSE)
colors.3<-setMap(Obj.Standard.Relaxed.Arboreal.simmap.ARD,c("gray84","lightskyblue"))
plot(colors.3,Obj.Standard.Relaxed.Arboreal.simmap.ARD,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")
title(main = "Standard.Relaxed.Arboreal.simmap.ARD")
#Gekkota clades
arc.cladelabels(text="Gekkota",node=findMRCA(StandardTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
#Dactyloidae
arc.cladelabels(text="Dactyloidae",node=findMRCA(StandardTree,c("Anolis_bonairensis","Anolis_trachyderma")),
                ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Scincidae
arc.cladelabels(text="Scincidae",node=3052, ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)


Obj.Standard.Relaxed.Arboreal.simmap.ER <-densityMap(Standard.Relaxed.Arboreal.simmap.ER,states=levels(Standard.Relaxed.Arboreal)[1:2],plot=FALSE)
colors.4<-setMap(Obj.Standard.Relaxed.Arboreal.simmap.ER,c("gray84","lightskyblue"))
plot(colors.4,Obj.Standard.Relaxed.Arboreal.simmap.ER,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Standard.Relaxed.Arboreal.simmap.ER")

## Label padbearing species on tips

#Standard.toepad.tips<-read.csv("Limbless_sse_dataset_strict_toepads.csv")
#Standard.toepad.tips <- filter(Standard.toepad.tips, toepads == "1")
#Standard.Toepad.Species <- Standard.toepad.tips$Concatenated
#ii<-sapply(Standard.Toepad.Species,grep,StandardTree$tip.label)
#ii
#ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")

#Examine the distribution of estimated origins and losses of arboreality in the Relaxed arboreality dataset

Standard.Relaxed.Arboreal.ARD.HPD<-describe.simmap(Standard.Relaxed.Arboreal.simmap.ARD)
Standard.Relaxed.Arboreal.ER.HPD<-describe.simmap(Standard.Relaxed.Arboreal.simmap.ER)

#Summary statistics

Standard.Relaxed.Arboreal.simmap.ARD.density <- density(Standard.Relaxed.Arboreal.simmap.ARD) #Summary statistics
Standard.Relaxed.Arboreal.simmap.ER.density <- density(Standard.Relaxed.Arboreal.simmap.ER) #Summary statistics


## ARD

#Relaxed arboreality frequency distribution
Standard.Relaxed.Arbor.ARD.origins.losses<- Standard.Relaxed.Arboreal.ARD.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Standard.Relaxed.Arbor.ARD.origins.losses) <- c("Gains", "Losses") #rename columns
Standard.Relaxed.Arbor.ARD.origins.losses.df <- as.data.frame(Standard.Relaxed.Arbor.ARD.origins.losses)

#gains of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)

Gains.Distribution.Standard.Relaxed.Arbor.ARD <- ggplot(Standard.Relaxed.Arbor.ARD.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+ggtitle("Gains.Distribution.Standard.Relaxed.Arbor.ARD")+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 255, xend = 300, 
               y = 25, yend = 25,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))
#losses of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Standard.Relaxed.Arbor.ARD <- ggplot(Standard.Relaxed.Arbor.ARD.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+ggtitle("Losses.Distribution.Standard.Relaxed.Arbor.ARD")+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 206, xend = 261, 
               y = 25, yend = 25,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))

#Put both histograms together
Standard.Relaxed.Arbor.ARD.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Relaxed.Arbor.ARD, Losses.Distribution.Standard.Relaxed.Arbor.ARD,labels = c("A", "B", ncol = 1, nrow = 1))

## ER

#Relaxed arboreality frequency distribution
Standard.Relaxed.Arbor.ER.origins.losses<- Standard.Relaxed.Arboreal.ER.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Standard.Relaxed.Arbor.ER.origins.losses) <- c("Gains", "Losses") #rename columns
Standard.Relaxed.Arbor.ER.origins.losses.df <- as.data.frame(Standard.Relaxed.Arbor.ER.origins.losses)
#gains of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Standard.Relaxed.Arbor.ER <- ggplot(Standard.Relaxed.Arbor.ER.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="lightskyblue", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+ggtitle("Gains.Distribution.Standard.Relaxed.Arbor.ER")+
  xlab("Arboreality gains") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 261, xend = 315, 
               y = 20, yend = 20,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))
#losses of Relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Standard.Relaxed.Arbor.ER <- ggplot(Standard.Relaxed.Arbor.ER.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="gray84", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+ggtitle("Losses.Distribution.Standard.Relaxed.Arbor.ER")+
  xlab("Arboreality losses") + ylab("Relative frequency")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 188, xend = 244, 
               y = 21, yend = 21,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 330))

#Put both histograms together
Standard.Relaxed.Arbor.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Relaxed.Arbor.ER, Losses.Distribution.Standard.Relaxed.Arbor.ER,labels = c("A", "B", ncol = 1, nrow = 1))

Standard.Arbor.ARD.ER.origins.losses.joined <- grid.arrange(arrangeGrob(Gains.Distribution.Standard.Strict.Arbor.ARD, 
                         Losses.Distribution.Standard.Strict.Arbor.ARD, 
                         Gains.Distribution.Standard.Strict.Arbor.ER, 
                         Losses.Distribution.Standard.Strict.Arbor.ER,
                         Gains.Distribution.Standard.Relaxed.Arbor.ARD, 
                         Losses.Distribution.Standard.Relaxed.Arbor.ARD, 
                         Gains.Distribution.Standard.Relaxed.Arbor.ER, 
                         Losses.Distribution.Standard.Relaxed.Arbor.ER, ncol = 2), # Second row with 2 plots in 2 different columns
             nrow = 1)  

##########################################################################
######                                                              ######
######   	                  Standard Multistate                     ######
######                                                              ######
##########################################################################

Standard.Multi.Dataset <- read.csv("2360_StandardEcologyLimbed_Multi.csv",row.names=1)

#Subset some characters for subsequent analyses

Standard.Multi.Arboreal <- as.factor(setNames(Standard.Multi.Dataset[,6],rownames(Standard.Multi.Dataset))) #toepads Standard dataset

#Assign prior probability on root node for non-arboreality

Multi.pi<-setNames(c(0,1,0),c("arboreal","not.arboreal", "semi.arboreal"))

#fit two Mk models

#fitARD.Multi.arbor.Standard<-fitMk(StandardTree,Standard.Multi.Arboreal,pi=Multi.pi,model="ARD")
#AIC(fitARD.Multi.arbor.Standard)
#
#fitER.Multi.arbor.Standard<<-fitMk(StandardTree,Standard.Multi.Arboreal,pi=Multi.pi,model="ER")
#AIC(fitER.Multi.arbor.Standard)
#
##Compare AIC values and determine best model
#
#aic.Multi.arbor.Standard<-setNames(sapply(list(fitARD.Multi.arbor.Standard,fitER.Multi.arbor.Standard),AIC), c("fitARD","fitER"))
#aic.Multi.arbor.Standard
#aic.w(aic.Multi.arbor.Standard)

## Standard Multi Arboreal reconstructions

Standard.Multi.Arboreal.simmap.ARD<-make.simmap(StandardTree,Standard.Multi.Arboreal,nsim=500,pi=Multi.pi,model="ARD")
Standard.Multi.Arboreal.simmap.ER<-make.simmap(StandardTree,Standard.Multi.Arboreal,nsim=500,pi=Multi.pi,model="ER")
Standard.Multi.Arboreal.simmap.SYM<-make.simmap(StandardTree,Standard.Multi.Arboreal,nsim=500,pi=Multi.pi,model="SYM")

#plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Standard.Multi.Arboreal.simmap.ARD<-summary(Standard.Multi.Arboreal.simmap.ARD)
cols<-setNames(c("#D8B365", "#F5F5F5", "#5AB4AC"),levels(Standard.Multi.Arboreal))
plot(Obj.Standard.Multi.Arboreal.simmap.ARD,fsize=0.6,ftype="off",colors=cols,ylim=c(-2,Ntip(StandardTree)))
add.simmap.legend(colors=cols[1:3],prompt=FALSE,x=0,y=2200,vertical=TRUE)
title(main = "Standard Multi Arbor ARD 500")
cladelabels(tree=StandardTree, "Gekkota", node=2362, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=StandardTree, "Dactyloidae", node= 4526, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=StandardTree, "Scincidae", node=3052, wing.length=NULL, cex=1,offset=1.5)
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey", offset = 1.4)

Obj.Standard.Multi.Arboreal.simmap.ER<-summary(Standard.Multi.Arboreal.simmap.ER)
cols<-setNames(c("#D8B365", "#F5F5F5", "#5AB4AC"),levels(Standard.Multi.Arboreal))
plot(Obj.Standard.Multi.Arboreal.simmap.ER,fsize=0.6,ftype="off",colors=cols,ylim=c(-2,Ntip(StandardTree)))
add.simmap.legend(colors=cols[1:3],prompt=FALSE,x=0,y=2200,vertical=TRUE)
title(main = "Standard Multi Arbor ER 500")
cladelabels(tree=StandardTree, "Gekkota", node=2362, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=StandardTree, "Dactyloidae", node= 4526, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=StandardTree, "Scincidae", node=3052, wing.length=NULL, cex=1,offset=1.5)
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey", offset = 1.4)

Obj.Standard.Multi.Arboreal.simmap.SYM<-summary(Standard.Multi.Arboreal.simmap.SYM)
cols<-setNames(c("#D8B365", "#F5F5F5", "#5AB4AC"),levels(Standard.Multi.Arboreal))
plot(Obj.Standard.Multi.Arboreal.simmap.SYM,fsize=0.6,ftype="off",colors=cols,ylim=c(-2,Ntip(StandardTree)))
add.simmap.legend(colors=cols[1:3],prompt=FALSE,x=0,y=2200,vertical=TRUE)
title(main = "Standard Multi Arbor SYM 500")
cladelabels(tree=StandardTree, "Gekkota", node=2362, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=StandardTree, "Dactyloidae", node= 4526, wing.length=NULL, cex=1,offset=1.5)
cladelabels(tree=StandardTree, "Scincidae", node=3052, wing.length=NULL, cex=1,offset=1.5)
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey", offset = 1.4)

#Examine the distribution of estimated origins and losses of arboreality in the Multi arboreality dataset

Standard.Multi.Arboreal.ARD.HPD<-describe.simmap(Standard.Multi.Arboreal.simmap.ARD)
Standard.Multi.Arboreal.ER.HPD<-describe.simmap(Standard.Multi.Arboreal.simmap.ER)
Standard.Multi.Arboreal.SYM.HPD<-describe.simmap(Standard.Multi.Arboreal.simmap.SYM)

#Summary statistics
Standard.Multi.Arboreal.simmap.ARD.density <- density(Standard.Multi.Arboreal.simmap.ARD) #Summary statistics
Standard.Multi.Arboreal.simmap.ER.density <- density(Standard.Multi.Arboreal.simmap.ER) #Summary statistics
Standard.Multi.Arboreal.simmap.SYM.density <- density(Standard.Multi.Arboreal.simmap.SYM) #Summary statistics

#Multi arboreality frequency distribution
## SYM multi distrib full
Standard.Multi.Arbor.SYM.origins.losses <- Standard.Multi.Arboreal.SYM.HPD$count[,2:7] #all 500 reconstructions with number of origins and reconstructions for each
colnames(Standard.Multi.Arbor.SYM.origins.losses) <- c("Arboreal.to.Non.arboreal", "Arboreal.to.Semi.arboreal","Non.arboreal.to.Arboreal","Non.arboreal.to.Semi.arboreal","Semi.arboreal.to.Arboreal","Semi.arboreal.to.Non.arboreal") #rename columns
Standard.Multi.Arbor.SYM.origins.losses <- as.data.frame(Standard.Multi.Arbor.SYM.origins.losses)
summary(Standard.Multi.Arbor.SYM.origins.losses) #Summary statistics for reconstructions
#HPD Calcs
Arboreal.to.Semi.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Standard.Multi.Arbor.SYM.origins.losses))[2,]
#110 - 156
Non.arboreal.to.Semi.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Standard.Multi.Arbor.SYM.origins.losses))[4,]
#298 - 359
Semi.arboreal.to.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Standard.Multi.Arbor.SYM.origins.losses))[5,]
#185 - 225
Semi.arboreal.to.Non.arboreal.HPD.SYM <- HPDinterval(as.mcmc(Standard.Multi.Arbor.SYM.origins.losses))[6,]
#167 - 222

#Arboreal to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.SYM.HPD.Arboreal.to.Semi.Arboreal <- ggplot(Standard.Multi.Arbor.SYM.origins.losses, aes(x=Arboreal.to.Semi.arboreal))+
  geom_histogram(color="#A6611A", fill="#A6611A",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Semi.arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.SYM")+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 110, xend = 156, 
               y = 27, yend = 27,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))
#Non-arboreality to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.SYM.HPD.NonArboreal.to.Semi.Arboreal <- ggplot(Standard.Multi.Arbor.SYM.origins.losses, aes(x=Non.arboreal.to.Semi.arboreal))+
  geom_histogram(color="#DFC27D", size=0.2,fill="#DFC27D",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Semi.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.SYM")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 298, xend = 359, 
               y = 17, yend = 17,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))


#Semi-arboreality to arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.SYM.HPD.SemiArboreal.to.Arboreal <- ggplot(Standard.Multi.Arbor.SYM.origins.losses, aes(x=Semi.arboreal.to.Arboreal))+
  geom_histogram(color="gray", size=0.2,fill="gray",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.SYM")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 185, xend = 225, 
               y = 28, yend = 28,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))

#Semi-arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.SYM.HPD.SemiArboreal.to.Non.Arboreal <- ggplot(Standard.Multi.Arbor.SYM.origins.losses, aes(x=Semi.arboreal.to.Non.arboreal))+
  geom_histogram(color="#80CDC1", size=0.2,fill="#80CDC1",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.SYM")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 167, xend = 222, 
               y = 25, yend = 25,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))


#Put everything together
Multistate.SYM.Standard <- grid.arrange(Standard.Multi.SYM.HPD.SemiArboreal.to.Non.Arboreal,
                                    Standard.Multi.SYM.HPD.SemiArboreal.to.Arboreal,
                                    Standard.Multi.SYM.HPD.NonArboreal.to.Semi.Arboreal,
                                    Standard.Multi.SYM.HPD.Arboreal.to.Semi.Arboreal, layout_matrix = rbind(c(1,1,1),c(2,3,4)))

Multistate.SYM.Standard <- ggarrange(Standard.Multi.SYM.HPD.SemiArboreal.to.Non.Arboreal,
                                     Standard.Multi.SYM.HPD.SemiArboreal.to.Arboreal,
                                     Standard.Multi.SYM.HPD.NonArboreal.to.Semi.Arboreal,
                                     Standard.Multi.SYM.HPD.Arboreal.to.Semi.Arboreal,labels = c("A", "B", "C", "D", ncol = 3, nrow = 3))

## ER multi distrib Standard
Standard.Multi.Arbor.ER.origins.losses <- Standard.Multi.Arboreal.ER.HPD$count[,2:7] #all 500 reconstructions with number of origins and reconstructions for each
colnames(Standard.Multi.Arbor.ER.origins.losses) <- c("Arboreal.to.Non.arboreal", "Arboreal.to.Semi.arboreal","Non.arboreal.to.Arboreal","Non.arboreal.to.Semi.arboreal","Semi.arboreal.to.Arboreal","Semi.arboreal.to.Non.arboreal") #rename columns
Standard.Multi.Arbor.ER.origins.losses <- as.data.frame(Standard.Multi.Arbor.ER.origins.losses)
summary(Standard.Multi.Arbor.ER.origins.losses) #Summary statistics for reconstructions
#HPD Calcs
Arboreal.to.Semi.arboreal.HPD.ER <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ER.origins.losses))[2,]
#84 - 115
Non.arboreal.to.Semi.arboreal.HPD.ER <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ER.origins.losses))[4,]
#239 - 281
Semi.arboreal.to.arboreal.HPD.ER <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ER.origins.losses))[5,]
#66 - 98
Semi.arboreal.to.Non.arboreal.HPD.ER <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ER.origins.losses))[6,]
#74 - 117
Arboreal.to.Non.arboreal.HPD.ER <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ER.origins.losses))[1,]
#48-76
Non.arboreal.to.Arboreal.HPD.ER <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ER.origins.losses))[3,]
#87-119
#Arboreal to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ER.HPD.Arboreal.to.Semi.Arboreal <- ggplot(Standard.Multi.Arbor.ER.origins.losses, aes(x=Arboreal.to.Semi.arboreal))+
  geom_histogram(color="#8C510A", fill="#8C510A",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Semi.arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ER")+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 84, xend = 115, 
               y = 32, yend = 32,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))
#Non-arboreality to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ER.HPD.NonArboreal.to.Semi.Arboreal <- ggplot(Standard.Multi.Arbor.ER.origins.losses, aes(x=Non.arboreal.to.Semi.arboreal))+
  geom_histogram(color="#D8B365", size=0.2,fill="#D8B365",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Semi.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 239, xend = 281, 
               y = 25, yend = 25,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))


#Semi-arboreality to arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ER.HPD.SemiArboreal.to.Arboreal <- ggplot(Standard.Multi.Arbor.ER.origins.losses, aes(x=Semi.arboreal.to.Arboreal))+
  geom_histogram(color="gray", size=0.2,fill="gray",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 66, xend = 98, 
               y = 35, yend = 35,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))

#Semi-arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ER.HPD.SemiArboreal.to.Non.Arboreal <- ggplot(Standard.Multi.Arbor.ER.origins.losses, aes(x=Semi.arboreal.to.Non.arboreal))+
  geom_histogram(color="#C7EAE5", size=0.2,fill="#C7EAE5",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 74, xend = 117, 
               y = 26, yend = 26,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))

#Arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ER.HPD.Arboreal.to.Non.Arboreal <- ggplot(Standard.Multi.Arbor.ER.origins.losses, aes(x=Arboreal.to.Non.arboreal))+
  geom_histogram(color="#5AB4AC", size=0.2,fill="#5AB4AC",binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Non.arboreal)),
                                                                                                             color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 48, xend = 76, 
               y = 37, yend = 37,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))

#Non-arboreality to Arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ER.HPD.Non.Arboreal.to.Arboreal <- ggplot(Standard.Multi.Arbor.ER.origins.losses, aes(x=Non.arboreal.to.Arboreal))+
  geom_histogram(color="#01665E", size=0.2,fill="#01665E",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Arboreal)),
                                                                                                             color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ER")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 87, xend = 119, 
               y = 30, yend = 30,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+coord_cartesian(xlim = c(0, 300))+ theme(plot.title = element_text(size=10))


#Put everything together

Multistate.ER.Standard <- ggarrange(Standard.Multi.ER.HPD.SemiArboreal.to.Non.Arboreal,
                                Standard.Multi.ER.HPD.Arboreal.to.Non.Arboreal, 
                                Standard.Multi.ER.HPD.SemiArboreal.to.Arboreal,
                                Standard.Multi.ER.HPD.NonArboreal.to.Semi.Arboreal,
                                Standard.Multi.ER.HPD.Arboreal.to.Semi.Arboreal, 
                                Standard.Multi.ER.HPD.Non.Arboreal.to.Arboreal,labels = c("A", "B", "C", "D", "E", "F", ncol = 3, nrow = 3))



## ARD multi distrib Standard
Standard.Multi.Arbor.ARD.origins.losses <- Standard.Multi.Arboreal.ARD.HPD$count[,2:7] #all 500 reconstructions with number of origins and reconstructions for each
colnames(Standard.Multi.Arbor.ARD.origins.losses) <- c("Arboreal.to.Non.arboreal", "Arboreal.to.Semi.arboreal","Non.arboreal.to.Arboreal","Non.arboreal.to.Semi.arboreal","Semi.arboreal.to.Arboreal","Semi.arboreal.to.Non.arboreal") #rename columns
Standard.Multi.Arbor.ARD.origins.losses <- as.data.frame(Standard.Multi.Arbor.ARD.origins.losses)
summary(Standard.Multi.Arbor.ARD.origins.losses) #Summary statistics for reconstructions
#HPD Calcs
Arboreal.to.Semi.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ARD.origins.losses))[2,]
#157 - 211
Non.arboreal.to.Semi.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ARD.origins.losses))[4,]
#223 - 280
Semi.arboreal.to.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ARD.origins.losses))[5,]
#154 - 193
Semi.arboreal.to.Non.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ARD.origins.losses))[6,]
#326 - 397
Arboreal.to.Non.arboreal.HPD.ARD <- HPDinterval(as.mcmc(Standard.Multi.Arbor.ARD.origins.losses))[1,]
#0-3
#Arboreal to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ARD.HPD.Arboreal.to.Semi.Arboreal <- ggplot(Standard.Multi.Arbor.ARD.origins.losses, aes(x=Arboreal.to.Semi.arboreal))+
  geom_histogram(color="#A6611A", fill="#A6611A",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Semi.arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ARD")+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 141, xend = 189, 
               y = 27, yend = 27,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))
#Non-arboreality to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ARD.HPD.NonArboreal.to.Semi.Arboreal <- ggplot(Standard.Multi.Arbor.ARD.origins.losses, aes(x=Non.arboreal.to.Semi.arboreal))+
  geom_histogram(color="#DFC27D", size=0.2,fill="#DFC27D",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Semi.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 261, xend = 318, 
               y = 24.5, yend = 24.5,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))


#Semi-arboreality to arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ARD.HPD.SemiArboreal.to.Arboreal <- ggplot(Standard.Multi.Arbor.ARD.origins.losses, aes(x=Semi.arboreal.to.Arboreal))+
  geom_histogram(color="gray", size=0.2,fill="gray",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 161, xend = 197, 
               y = 27, yend = 27,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))

#Semi-arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ARD.HPD.SemiArboreal.to.Non.Arboreal <- ggplot(Standard.Multi.Arbor.ARD.origins.losses, aes(x=Semi.arboreal.to.Non.arboreal))+
  geom_histogram(color="#80CDC1", size=0.2,fill="#80CDC1",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 296, xend = 367, 
               y = 17, yend = 17,
               colour = "black")+scale_y_continuous(limits = c(0,30), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))

#Arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
Standard.Multi.ARD.HPD.Arboreal.to.Non.Arboreal <- ggplot(Standard.Multi.Arbor.ARD.origins.losses, aes(x=Arboreal.to.Non.arboreal))+
  geom_histogram(color="#018571", size=0.2,fill="#018571",binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Non.arboreal)),
                                                                                                             color="black", linetype="dashed", size=1)+ggtitle("Multi.Distribution.Standard.Arbor.ARD")+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 1, xend = 7, 
               y = 125, yend = 125,
               colour = "black")+scale_y_continuous(limits = c(0,200), expand = c(0, 0))+coord_cartesian(xlim = c(0, 400))+ theme(plot.title = element_text(size=10))



#Put everything together
Multistate.ARD.Standard <- grid.arrange(Standard.Multi.ARD.HPD.Arboreal.to.Non.Arboreal,
                                    Standard.Multi.ARD.HPD.SemiArboreal.to.Non.Arboreal,
                                    Standard.Multi.ARD.HPD.SemiArboreal.to.Arboreal,
                                    Standard.Multi.ARD.HPD.NonArboreal.to.Semi.Arboreal,
                                    Standard.Multi.ARD.HPD.Arboreal.to.Semi.Arboreal, layout_matrix = rbind(c(1,1,1,1),c(2,3,4,5)))



##########################################################################
######                                                              ######
######   	                  Standard Toepads                        ######
######                                                              ######
##########################################################################

Standard.Toepads <- as.factor(setNames(Standard.Strict.Dataset[,2],rownames(Standard.Strict.Dataset))) #toepads Standard dataset

#Assign prior probability on root node for non-arboreality

toepads.pi<-setNames(c(1,0),c("0","1"))

#fit two Mk models

#fitARD.Toepads.Standard<-fitMk(StandardTree,Standard.Toepads,pi=toepads.pi,model="ARD")
#AIC(fitARD.Toepads.Standard)
#
#fitER.Toepads.Standard<<-fitMk(StandardTree,Standard.Toepads,pi=toepads.pi,model="ER")
#AIC(fitER.Toepads.Standard)
#
##Compare AIC values and determine best model
#
#aic.Toepads.Standard<-setNames(sapply(list(fitARD.Toepads.Standard,fitER.Toepads.Standard),AIC), c("fitARD","fitER"))
#aic.Toepads.Standard
#aic.w(aic.Toepads.Standard)

## Standard Strict Arboreal reconstructions

Standard.Toepads.simmap.ARD<-make.simmap(StandardTree,Standard.Toepads,nsim=500,pi=toepads.pi,model="ARD")
Standard.Toepads.simmap.ER<-make.simmap(StandardTree,Standard.Toepads,nsim=500,pi=toepads.pi,model="ER")

#plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

Obj.Standard.Toepads.simmap.ARD <-densityMap(Standard.Toepads.simmap.ARD,states=levels(Standard.Toepads)[1:2],plot=FALSE)
colors.5<-setMap(Obj.Standard.Toepads.simmap.ARD,c("gray84","lightskyblue"))
plot(colors.5,Obj.Standard.Toepads.simmap.ARD,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Standard Toepads ARD 500")
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")
arc.cladelabels(text="Gekkota",node=findMRCA(StandardTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
arc.cladelabels(text="Dactyloidae",node=findMRCA(StandardTree,c("Anolis_bonairensis","Anolis_trachyderma")),ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
arc.cladelabels(text="Scincidae",node=3052, ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)


Obj.Standard.Toepads.simmap.ER <-densityMap(Standard.Toepads.simmap.ER,states=levels(Standard.Toepads)[1:2],plot=FALSE)
colors.6<-setMap(Obj.Standard.Toepads.simmap.ER,c("gray84","lightskyblue"))
plot(colors.6,Obj.Standard.Toepads.simmap.ER,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)
title(main = "Standard Toepads ER 500")
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")
arc.cladelabels(text="Gekkota",node=findMRCA(StandardTree,c("Phyllurus_kabikabi","Phelsuma_lineata")),ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
arc.cladelabels(text="Dactyloidae",node=findMRCA(StandardTree,c("Anolis_bonairensis","Anolis_trachyderma")),ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)
arc.cladelabels(text="Scincidae",node=3052, ln.offset=1.03,lab.offset=1.06, mark.node=FALSE)

#Label strict arboreal species on tips
#Standard.Strict.Arboreal.tips<-read.csv("Limbless_sse_dataset_strict_toepads.csv")
#Standard.Strict.Arboreal.tips <- filter(Standard.Strict.Arboreal.tips, strict.arboreal == "1")
#Standard.Strict.Arboreal.tips.Species <- Standard.Strict.Arboreal.tips$Concatenated
#arbor.ii<-sapply(Standard.Strict.Arboreal.tips.Species,grep,StandardTree$tip.label)
#arbor.ii
#ape::tiplabels(pch=21, tip=arbor.ii, cex=2, col='black',bg="grey")

#Examine the distribution of estimated origins and losses of arboreality in the strict arboreality dataset

Standard.Toepads.ARD.HPD<-describe.simmap(Standard.Toepads.simmap.ARD)
Standard.Toepads.ER.HPD<-describe.simmap(Standard.Toepads.simmap.ER)

#Summary statistics
Standard.Toepads.simmap.ARD.density <- density(Standard.Toepads.simmap.ARD) #Summary statistics
Standard.Toepads.simmap.ER.density <- density(Standard.Toepads.simmap.ER) #Summary statistics

##Strict arboreality frequency distribution
Standard.Toepads.ARD.origins.losses<- Standard.Toepads.ARD.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Standard.Toepads.ARD.origins.losses) <- c("Gains", "Losses") #rename columns
Standard.Toepads.ARD.origins.losses.df <- as.data.frame(Standard.Toepads.ARD.origins.losses)
#gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Standard.Toepads.ARD <- ggplot(Standard.Toepads.ARD.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+ggtitle("Gains.Distribution.Standard.Toepads.ARD")+
  xlab("Toepad gains") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 3, xend = 5, 
               y = 350, yend = 350,
               colour = "black")+scale_y_continuous(limits = c(0,360), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30))
#losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Standard.Toepads.ARD <- ggplot(Standard.Toepads.ARD.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="black", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+ggtitle("Losses.Distribution.Standard.Toepads.ARD")+
  xlab("Toepad losses") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 20, xend = 26, 
               y = 200, yend = 200,
               colour = "black")+scale_y_continuous(limits = c(0,360), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30))

#Put both histograms together
Standard.Toepads.ARD.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Toepads.ARD, Losses.Distribution.Standard.Toepads.ARD,labels = c("A", "B", ncol = 1, nrow = 1))

#Strict arboreality frequency distribution
Standard.Toepads.ER.origins.losses<- Standard.Toepads.ER.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(Standard.Toepads.ER.origins.losses) <- c("Gains", "Losses") #rename columns
Standard.Toepads.ER.origins.losses.df <- as.data.frame(Standard.Toepads.ER.origins.losses)
#gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Gains.Distribution.Standard.Toepads.ER <- ggplot(Standard.Toepads.ER.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+ggtitle("Gains.Distribution.Standard.Toepads.ER")+
  xlab("Toepad gains") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 8, xend = 11, 
               y = 325, yend = 325,
               colour = "black")+scale_y_continuous(limits = c(0,350), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30))
#losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
Losses.Distribution.Standard.Toepads.ER <- ggplot(Standard.Toepads.ER.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="black", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+ggtitle("Losses.Distribution.Standard.Toepads.ER")+
  xlab("Toepad losses") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 11, xend = 14, 
               y = 260, yend = 260,
               colour = "black")+scale_y_continuous(limits = c(0,350), expand = c(0, 0))+coord_cartesian(xlim = c(0, 30))

#Put both histograms together
Standard.Toepads.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Toepads.ER, Losses.Distribution.Standard.Toepads.ER,labels = c("A", "B", ncol = 1, nrow = 1))

##All

Standard.Toepads.ARD.ER.origins.losses.joined <- ggarrange(Gains.Distribution.Standard.Toepads.ARD, Losses.Distribution.Standard.Toepads.ARD,Gains.Distribution.Standard.Toepads.ER,Losses.Distribution.Standard.Toepads.ER,labels = c("A", "B", "C", "D",ncol = 2, nrow = 2))


