
##################################################################################################
######                                                                                      ######
###### "Novel Tests of the Key Innovation Hypothesis: Adhesive Toepads in Arboreal Lizards" ######
######                      Authors: Aryeh H. Miller and James T. Stroud                    ######
######                                                                                      ######
##################################################################################################
  
#Load relevant packages

library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(coda)
library(tidyverse)
library(plyr)
library(ggpubr)
library(scales)
library(tidybayes)
library(dplyr)
library(plotrix)
library(cowplot)
library(RColorBrewer)

#Read time-calibrated phylogeny of 2376 lizard species-- see Ramm et al. (2020) for details regarding construction and calibration methods.

tree <- read.tree(file = "2376Lizards.tre")

#Read in ecological and morphological (toepad) trait dataset (Supplementary Table 1) 
#Name rows with concatenated genus_species column

data<-read.csv("SupplementaryTable1_EcologicalTraitDatabase_September2020.csv",row.names=3)

#Subset some characters for subsequent analyses

strict_arboreal<-as.factor(setNames(data[,3],rownames(data))) #Strict arboreal dataset
multi_arboreal<-as.factor(setNames(data[,4],rownames(data))) #Multistate dataset
relaxed_arboreal<-as.factor(setNames(data[,5],rownames(data))) #Relaxed arboreal dataset
toepads<-as.factor(setNames(data[,6],rownames(data))) #toepad dataset

#NOTE: if you don't wait for the computation with the analyses, simply load the relevant .RData files. These analyses, specifically the SCM and the evolutionary correlation analyses will take a bit.

##########################################################################
######                                                              ######
######   	Strict arboreality stochastic character mapping           ######
######                                                              ######
##########################################################################

#Assign prior probability on root node for non-arboreality

pi<-setNames(c(0,1),c("arboreal","not.arboreal"))
pi

#fit two Mk models

fitARD<-fitMk(tree,strict_arboreal,pi=pi,model="ARD")
AIC(fitARD)

fitER<-fitMk(tree,strict_arboreal,pi=pi,model="ER")
AIC(fitER)

#Compare AIC values and determine best model

aic<-setNames(sapply(list(fitARD,fitER),AIC), c("fitARD","fitER"))
aic
aic.w(aic)

#Use ARD

strict_arboreal_simmap<-make.simmap(tree,strict_arboreal,nsim=500,pi=pi,model="ARD")

#plot posterior distribution of changes on the tree-- set Carolina blue and grey color scheme.

obj_strict_arboreal<-densityMap(strict_arboreal_simmap,states=levels(strict_arboreal)[2:1],plot=FALSE)
colors<-setMap(obj_strict_arboreal,c("gray84","lightskyblue"))
plot(colors,obj_strict_arboreal,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(2,7),legend=79.535)

#Label padbearing species on tips
Toepad.tips<-read.csv("SupplementaryTable1_EcologicalTraitDatabase_September2020.csv")
Toepad.tips <- filter(Toepad.tips, toepads == "1")
Toepad.species <- Toepad.tips$Concatenated
ii<-sapply(Toepad.species,grep,tree$tip.label)
ii
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")

#Examine the distribution of estimated origins and losses of arboreality in the strict arboreality dataset

strict.arboreal.HPD<-describe.simmap(strict_arboreal_simmap)
strict.arboreal.HPD
#Summary statistics
density.strict.arboreal.empirical <- density(strict_arboreal_simmap) #Summary statistics

#Strict arboreality frequency distribution
strictarbor.origins.losses<- strict.arboreal.HPD$count[,2:3] #all 500 reconstructions with number of origins and losses for each
colnames(strictarbor.origins.losses) <- c("Losses", "Gains") #rename columns
strictarbor_origins_losses_df <- as.data.frame(strictarbor.origins.losses)
#gains of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
gains_distribution_strictarboreality <- ggplot(strictarbor_origins_losses_df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality gains") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 94, xend = 122, 
               y = 37, yend = 37,
               colour = "black")+scale_y_continuous(limits = c(0,40), expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 261, by = 5))
#losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
losses_distribution_strictarboreality <- ggplot(strictarbor_origins_losses_df, aes(x=Losses))+
  geom_histogram(color="black", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+
  xlab("Arboreality losses") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 190, xend = 250, 
               y = 19.5, yend = 19.5,
               colour = "black")+scale_y_continuous(limits = c(0,20.5), expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 265, by = 10))

#Put both histograms together
StrictArboreality_OriginsLosses_Frequency_Plot <- ggarrange(gains_distribution_strictarboreality, losses_distribution_strictarboreality,labels = c("A", "B", ncol = 1, nrow = 1))


##########################################################################
######                                                              ######
######   	Relaxed arboreality stochastic character mapping          ######
######                                                              ######
##########################################################################

### SCM Arboreal RELAXED
#Set prior on root node for non-arboreality
pi.arboreal.relaxed<-setNames(c(0,1),c("arboreal","not.arboreal"))
pi.arboreal.relaxed

#fitmk relaxed
fitARD.relaxed<-fitMk(tree,relaxed_arboreal,pi=pi.arboreal.relaxed,model="ARD")
AIC(fitARD.relaxed)
fitER.relaxed<-fitMk(tree,relaxed_arboreal,pi=pi.arboreal.relaxed,model="ER")
AIC(fitER.relaxed)
#aicw
aic<-setNames(sapply(list(fitARD.relaxed,fitER.relaxed),AIC), c("fitARD","fitER"))
aic
aic.w(aic)

relaxed.arboreality.SCM <-make.simmap(tree,relaxed_arboreal,nsim=500,pi=pi.arboreal.relaxed,model="ARD")
summary.relaxed.arboreality<-summary(relaxed.arboreality.SCM)
summary.relaxed.arboreality

#posterior distribution of changes on the tree (rather than states)
obj.relaxed.arboreality<-densityMap(relaxed.arboreality.SCM,states=levels(relaxed_arboreal)[2:1],plot=FALSE)
colors<-setMap(obj.relaxed.arboreality,c("gray84","lightskyblue"))
plot(colors,obj.relaxed.arboreality,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1, lwd=c(1,7),legend=79.535)

#Label padbearing species on tips
ape::tiplabels(pch=21, tip=ii, cex=2, col='black',bg="grey")

#relaxed arboreality frequency distribution

density.relaxed.arboreal <- density(relaxed.arboreality.SCM)
relaxed.arboreal.HPD<-describe.simmap(relaxed.arboreality.SCM)
relaxed.arboreal.HPD

relaxedarbor.origins.losses<- relaxed.arboreal.HPD$count[,2:3] #all 500 reconstructions with number of origins and reconstructions for each
colnames(relaxedarbor.origins.losses) <- c("Losses", "Gains") #rename columns
relaxedarbor_origins_losses_df <- as.data.frame(relaxedarbor.origins.losses)
#gains of relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
gains_distribution_relaxedarboreality <- ggplot(relaxedarbor_origins_losses_df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality gains") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                                color="black", linetype="dashed", size=1)+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 228, xend = 273, 
               y = 25, yend = 25,
               colour = "black")+scale_y_continuous(limits = c(0,27), expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 295, by = 5))
#losses of relaxed arboreality with specified 95% HPD (solid line) and mean (dashed line)
losses_distribution_relaxedrboreality <- ggplot(relaxedarbor_origins_losses_df, aes(x=Losses))+
  geom_histogram(color="black", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+
  xlab("Arboreality losses") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                 color="black", linetype="dashed", size=1)+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 182, xend = 234, 
               y = 21.5, yend = 21.5,
               colour = "black")+scale_y_continuous(limits = c(0,23), expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 265, by = 5))

#Put both histograms together
RelaxedArboreality_OriginsLosses_Frequency_Plot <- ggarrange(gains_distribution_relaxedarboreality, losses_distribution_relaxedrboreality,labels = c("A", "B", ncol = 1, nrow = 1))

##########################################################################
######                                                              ######
######   	  Multistate habitat stochastic character mapping         ######
######                                                              ######
##########################################################################


### SCM Multistate Arboreality
multi_arboreal<-as.factor(setNames(data[,4],rownames(data))) #Multistate dataset
states_arboreal_multi<-sort(unique(multi_arboreal))

#Set prior on root node for non-arboreality

pi.arboreal.multi<-setNames(c(0,1,0),c("arboreal","not.arboreal", "semi.arboreal"))
pi.arboreal.multi

#fitmk Multistate

fitARD.multi<-fitMk(tree,arboreal_multi,pi=pi.arboreal.multi,model="ARD")
AIC(fitARD.multi)

fitER.multi<-fitMk(tree,arboreal_multi,pi=pi.arboreal.multi,model="ER")
AIC(fitER.multi)

fitSYM.multi<-fitMk(tree,arboreal_multi,pi=pi.arboreal.multi,model="SYM")
AIC(fitSYM.multi)

#AIC and AICw comparisons
aic<-setNames(sapply(list(fitARD.multi,fitER.multi,fitSYM.multi),AIC), c("fitARD","fitER","fitSYM"))
aic
aic.w(aic)

## multi SCM
multi.arboreality.SCM <-make.simmap(tree,multi_arboreal,nsim=500,pi=pi.arboreal.multi,model="ARD")
summary.multi.arboreality<-summary(multi.arboreality.SCM)
#colors
show_col(viridis_pal()(3))
#Pies on nodes
cols<-setNames(c("#440154FF","#21908CFF","#FDE725FF"),levels(states_arboreal_multi))
par(fg="transparent",lend=1)
plot(summary.multi.arboreality,fsize=0,ftype="off",ln=0.4,type="fan",edge.width = 0.5,offset=1,colors=cols)
add.simmap.legend(colors=cols[3:1],prompt=FALSE,x=0,y=0,vertical=TRUE)
legend(x="bottomleft",levels(states_arboreal_multi),pch=22,
       pt.bg=cols,pt.cex=1.5,bty="n",cex=0.7)

#View probabilities at nodes if curious
View(summary.multi.arboreality$ace)

#pies
cols<-setNames(c("#440154FF","#21908CFF","#FDE725FF"),levels(states_arboreal_multi))
par(fg="transparent",lend=1)
plot(summary.multi.arboreality,fsize=0,ftype="off",ln=0.4,type="fan",edge.width = 0.5,offset=1,colors=cols)
add.simmap.legend(colors=cols[3:1],prompt=FALSE,x=0,y=0,vertical=TRUE)

#HPD Calculations and Frequency Distributions
multi.arboreal.HPD<-describe.simmap(multi.arboreality.SCM.1)
multi.arboreal.HPD
###strict arboreality frequency distribution
multiarbor.origins.losses<- multi.arboreal.HPD$count[,2:7] #all 500 reconstructions with number of origins and reconstructions for each
colnames(multiarbor.origins.losses) <- c("Arboreal.to.Non.arboreal", "Arboreal.to.Semi.arboreal","Non.arboreal.to.Arboreal","Non.arboreal.to.Semi.arboreal","Semi.arboreal.to.Arboreal","Semi.arboreal.to.Non.arboreal") #rename columns
multiarbor_origins_losses_df <- as.data.frame(multiarbor.origins.losses)
#HPD Calcs
Arboreal.to.Semi.arboreal.HPD <- HPDinterval(as.mcmc(multiarbor_origins_losses_df))[2,]
#176 - 244
Non.arboreal.to.Semi.arboreal.HPD <- HPDinterval(as.mcmc(multiarbor_origins_losses_df))[4,]
#205 - 256
Semi.arboreal.to.arboreal.HPD <- HPDinterval(as.mcmc(multiarbor_origins_losses_df))[5,]
#135 - 175
Semi.arboreal.to.Non.arboreal.HPD <- HPDinterval(as.mcmc(multiarbor_origins_losses_df))[6,]
#312 - 381

#Arboreal to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
HPD_distribution_Arboreal_to_Semi_arboreal <- ggplot(multiarbor_origins_losses_df, aes(x=Arboreal.to.Semi.arboreal))+
  geom_histogram(color="black", fill="#FDE725FF",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Arboreal.to.Semi.arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+
  geom_text(x = 250,  y = 26, 
            label = "", 
            colour = "black")+
  geom_segment(x = 176, xend = 244, 
               y = 19, yend = 19,
               colour = "black")+scale_y_continuous(limits = c(0,20), expand = c(0, 0))+ scale_x_continuous(limits = c(0, 410),expand= c(0,0))
#Non-arboreality to semi-arboreality with specified 95% HPD (solid line) and mean (dashed line)
HPD_distribution_NonArboreal_to_Semiarboreal <- ggplot(multiarbor_origins_losses_df, aes(x=Non.arboreal.to.Semi.arboreal))+
  geom_histogram(color="black", size=0.2,fill="#FDE725FF",binwidth=0.5)+theme_bw()+
  xlab("Non-arboreality to Semi-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Non.arboreal.to.Semi.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 205, xend = 256, 
               y = 23, yend = 23,
               colour = "black")+scale_y_continuous(limits = c(0,25), expand = c(0, 0))+ scale_x_continuous(limits = c(0, 410),expand= c(0,0))


#Semi-arboreality to arboreality with specified 95% HPD (solid line) and mean (dashed line)
HPD_distribution_SemiArboreal_to_Arboreal <- ggplot(multiarbor_origins_losses_df, aes(x=Semi.arboreal.to.Arboreal))+
  geom_histogram(color="black", size=0.2,fill="#440154FF",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Arboreal)),
                                                                                                              color="black", linetype="dashed", size=1)+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 135, xend = 175, 
               y = 32, yend = 32,
               colour = "black")+scale_y_continuous(limits = c(0,33), expand = c(0, 0))+ scale_x_continuous(limits = c(0, 410),expand= c(0,0))

#Semi-arboreality to Non-arboreality with specified 95% HPD (solid line) and mean (dashed line)
HPD_distribution_SemiArboreal_to_Non_arboreal <- ggplot(multiarbor_origins_losses_df, aes(x=Semi.arboreal.to.Non.arboreal))+
  geom_histogram(color="black", size=0.2,fill="#21908CFF",binwidth=0.5)+theme_bw()+
  xlab("Semi-arboreality to Non-arboreality") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Semi.arboreal.to.Non.arboreal)),
                                                                                                                  color="black", linetype="dashed", size=1)+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 312, xend = 381, 
               y = 18, yend = 18,
               colour = "black")+scale_y_continuous(limits = c(0,20), expand = c(0, 0))+ scale_x_continuous(limits = c(0, 410),expand= c(0,0))


#Put everything together
Multistate_Frequency_Plot <- ggarrange(HPD_distribution_Arboreal_to_Semi_arboreal, HPD_distribution_NonArboreal_to_Semiarboreal,HPD_distribution_SemiArboreal_to_Arboreal,HPD_distribution_SemiArboreal_to_Non_arboreal,labels = c("A", "B", "C", "D", ncol = 2, nrow = 2))

##########################################################################
######                                                              ######
######   	        Toepads stochastic character mapping              ######
######                                                              ######
##########################################################################


#Set prior on root node for no toepads
pi_toepads<-setNames(c(1,0),c("0","1"))
pi_toepads
#reconstruct toepad evolution
fitARD<-fitMk(tree,toepads,pi=pi_toepads,model="ARD")
AIC(fitARD)
fitER<-fitMk(tree,toepads,pi=pi_toepads,model="ER")
AIC(fitER)
#aicw
aic<-setNames(sapply(list(fitARD,fitER),AIC), c("fitARD","fitER"))
aic
aic.w(aic)

#SCM
toepads_500_simmap_ARD <-make.simmap(tree,toepads,nsim=500,pi=pi_toepads,model="ARD")
sum.toepads.updated<-summary(toepads_500_simmap_ARD)
Toepads.500.HPD <- describe.simmap(toepads_500_simmap_ARD)
toepads_500_simmap_ARD <- toepads_500_simmap_ARD_Sept2020_flipped

#posterior distribution of changes on the tree (rather than states)
obj_toepads <-densityMap(toepads_500_simmap_ARD,states=levels(toepads)[1:2],plot=FALSE)
colors<-setMap(obj_toepads,c("gray84","lightskyblue"))
plot(colors,obj_toepads,lwd=1,outline=F,fsize=c(0.7,0.9),type="fan",ftype="off", offset=1)

toepads.SCM.intervals.stats <-density(toepads_500_simmap_ARD) #summary stats
toepads.SCM.intervals.stats
pad.origins.losses<- Toepads.500.HPD$count[,2:3] #all 500 reconstructions with number of origins and reconstructions for each
colnames(pad.origins.losses) <- c("Gains", "Losses") #rename columns
pad_origins_losses_df <- as.data.frame(pad.origins.losses)
#gains with specified 95% HPD (solid line) and mean (dashed line)
gains_distribution <- ggplot(pad_origins_losses_df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",binwidth=0.5)+theme_bw()+
  xlab("Toepad gains") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Gains)),
                                                                                           color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 320, 
            label = "", 
            colour = "black")+
  geom_segment(x = 8, xend = 10, 
               y = 320, yend = 320,
               colour = "black")+scale_y_continuous(limits = c(0,350), expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 12, by = 1))
#losses with specified 95% HPD (solid line) and mean (dashed line)
losses_distribution <- ggplot(pad_origins_losses_df, aes(x=Losses))+
  geom_histogram(color="black", fill="gray84",binwidth=0.5)+theme_bw()+
  xlab("Toepad losses") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                            color="black", linetype="dashed", size=1)+
  geom_text(x = 12,  y = 209, 
            label = "", 
            colour = "black")+
  geom_segment(x = 9, xend = 14, 
               y = 202, yend = 202,
               colour = "black")+scale_y_continuous(limits = c(0,220), expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 20, by = 1))

#Put both histograms together
Toepads_OriginsLosses_Frequency_Plot <- ggarrange(gains_distribution, losses_distribution,labels = c("A", "B", ncol = 1, nrow = 1))


#Gekkota Toepad SCM Frequency Distributions
#Extract the 'Gekkota' clade
Gekkota_trees<-lapply(toepads_500_simmap_ARD_22June20,extract.clade.simmap,node=2385)
class(Gekkota_trees)<-"multiPhylo"
Gekkota_trees
Gekkota_summary <- describe.simmap(Gekkota_trees)
toepads_500_simmap_ARD_Gekkota<-densityMap(Gekkota_trees,states=levels(toepads)[1:2],plot=FALSE)
colors<-setMap(toepads_500_simmap_ARD_Gekkota,c("gray84","lightskyblue"))
plot(colors,toepads_500_simmap_ARD_Gekkota,outline=F,fsize=c(0.7,0.9),ftype="off", offset=1, lwd=c(1,7))

#Gekkota Toepad frequency distribution
Gekkota.toes.origins.losses<- Gekkota_summary$count[,2:3] #all 500 reconstructions with number of origins and reconstructions for each
colnames(Gekkota.toes.origins.losses) <- c("Gains", "Losses") #rename columns
Gekkota.toes.origins.losses.df <- as.data.frame(Gekkota.toes.origins.losses)
##
View(Gekkota.toes.origins.losses.df)
#95% HPD Interval
Gekkota.HPD.Gains <- HPDinterval(as.mcmc(Gekkota.toes.origins.losses.df))[1,]
Gekkota.HPD.Losses <- HPDinterval(as.mcmc(Gekkota.toes.origins.losses.df))[2,]

#gains of toepads with specified 95% HPD (solid line) and mean (dashed line)
gains_distribution_toepads_Gekkota<- ggplot(Gekkota.toes.origins.losses.df, aes(x=Gains))+
  geom_histogram(color="black", fill="lightskyblue",size=0.2,binwidth=0.5)+theme_bw()+
  xlab("Toepad gains in Gekkota") + ylab("Relative frequency across 500 stochastic maps")+
  geom_text(x = 1,  y = 500, 
            label = "", 
            colour = "black")+
  geom_segment(x = 6, xend = 8, 
               y = 375, yend = 375,
               colour = "black")+geom_vline(aes(xintercept=mean(Gains)),
                                            color="black", linetype="dashed", size=1)+scale_y_continuous(limits = c(0,400), expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 10, by = 1))
#losses of strict arboreality with specified 95% HPD (solid line) and mean (dashed line)
losses_distribution_toepads_Gekkota <- ggplot(Gekkota.toes.origins.losses.df, aes(x=Losses))+
  geom_histogram(color="black", size=0.2,fill="gray84",binwidth=0.5)+theme_bw()+
  xlab("Toepad losses in Gekkota") + ylab("Relative frequency across 500 stochastic maps")+ geom_vline(aes(xintercept=mean(Losses)),
                                                                                                       color="black", linetype="dashed", size=1)+
  geom_text(x = 200,  y = 22, 
            label = "", 
            colour = "black")+
  geom_segment(x = 8, xend = 12, 
               y = 225, yend = 225,
               colour = "black")+scale_y_continuous(limits = c(0,250), expand = c(0, 0))+ scale_x_continuous(breaks = seq(8, 16, by = 1))

#Put them together
Toepad_OriginsLosses_Frequency_Plot_Gekkota<- ggarrange(gains_distribution_toepads_Gekkota, losses_distribution_toepads_Gekkota,labels = c("A", "B", ncol = 2, nrow =2))

##########################################################################
######                                                              ######
######   	Stacked barplot multistate arboreality figure             ######
######                                                              ######
##########################################################################

#multistate without numbered proportions
data<-read.csv("SupplementaryTable1_EcologicalTraitDatabase_September2020.csv",row.names = 4)
#set up plot
family <- data$Family
arboreal <- data$strict.dataset
multi <- data$multistate
# execute plot (alphabetical order)
p<- ggplot(data, aes(family, fill=factor(multi,levels=c("arboreal","semi.arboreal","not.arboreal"))))+geom_bar(position="fill")+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+ scale_fill_manual(values=c("lightskyblue","gray71","grey89"))+scale_y_continuous(limits = c(0,1), expand = c(0, 0))+
  theme(panel.grid = element_blank(),panel.border = element_blank())+coord_flip()+labs(x ="Family", y = "Proportion of Species")+scale_x_discrete(limits = rev(levels(data$Family)))+ theme_dark()+ guides(fill=guide_legend(title="State"))

#Re-order given phylogenetic relationships from Pyron et al. 2013
family_phyloorder_list <- c("Dactyloidae",
                            "Corytophanidae",
                            "Liolaemidae",
                            "Leiosauridae",	
                            "Opluridae"	,
                            "Hoplocercidae"	,
                            "Polychrotidae",
                            "Phrynosomatidae",
                            "Crotaphytidae",
                            "Leiocephalidae",
                            "Iguanidae",
                            "Tropiduridae",
                            "Agamidae",	
                            "Chamaeleonidae",
                            "Varanidae",
                            "Lanthanotidae",
                            "Shinisauridae",
                            "Anguidae",
                            "Anniellidae",
                            "Helodermatidae",
                            "Xenosauridae",
                            "Lacertidae",
                            "Amphisbaenidae",
                            "Trogonophiidae",
                            "Cadeidae"	,
                            "Blanidae",
                            "Bipedidae",
                            "Rhineuridae",
                            "Gymnophthalmidae",
                            "Teiidae",
                            "Scincidae",
                            "Cordylidae",
                            "Gerrhosauridae",
                            "Xantusiidae",
                            "Gekkonidae",
                            "Phyllodactylidae",	
                            "Sphaerodactylidae",
                            "Eublepharidae",
                            "Diplodactylidae",
                            "Pygopodidae",
                            "Carphodactylidae",	
                            "Dibamidae" )
reordered_phylo <- rev(family_phyloorder_list)
p+xlim(reordered_phylo)

##########################################################################
######                                                              ######
######   	Evolutionary correlation and transition rate analyses     ######
######                                                              ######
##########################################################################

#Fitpagel with fixed ("fossilized") root state

#STRICT

#add a tip of zero-length to the root using "bind.tip" and then assigning it non-arboreal and padless
tree_root <-bind.tip(tree,tip.label="ROOT",edge.length=0,
                     where=Ntip(tree)+1)
plotTree(tree_root,ftype="off",lwd=1)
nn<-which(tree_root$tip.label=="ROOT")
tiplabels(tree_root$tip.label[nn],nn,adj=c(-0.1,0.5),
          frame="none",cex=0.8,font=3)
# add fossilized root states
arboreal_strict<-as.factor(setNames(data[,3],rownames(data)))
toepads<-as.factor(setNames(data[,6],rownames(data)))

arboreal_fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                  setNames(as.character(arboreal_strict),names(arboreal_strict))))
padless_fossilized <-as.factor(c(setNames("0","ROOT"),
                                  setNames(as.character(toepads),names(toepads))))

#fit fixed-root models for Strict dataset

fit.fixed_dep_y_ARD<-fitPagel(tree_root,padless_fossilized,arboreal_fossilized,dep.var="y",model="ARD") #strict
fit.fixed_dep_y_ARD
plot(fit.fixed_dep_y_ARD)
title("fit.fixed_dep_y_ARD_Strict")

fit.fixed_dep_x_ARD_Strict<-fitPagel(tree_root,padless_fossilized,arboreal_fossilized,dep.var="x",model="ARD") #strict
plot(fit.fixed_dep_x_ARD_Strict)
title("fit.fixed_dep_x_ARD_Strict")

fit.fixed_dep_xy_ARD_Strict<-fitPagel(tree_root,padless_fossilized,arboreal_fossilized,dep.var="xy",model="ARD") #strict
plot(fit.fixed_dep_xy_ARD_Strict)
title("fit.fixed_dep_xy_ARD_Strict")

aic_strict<-setNames(c(fit.fixed_dep_x_ARD_Strict$dependent.AIC,
                       fit.fixed_dep_y_ARD$dependent.AIC,
                       fit.fixed_dep_xy_ARD_Strict$dependent.AIC),
                     c("dependent x",
                       "dependent y","dependent x&y"))
aic_strict
aic.w(aic_strict)

#Relaxed fitpagel

arboreal_relaxed<-as.factor(setNames(data[,6],rownames(data)))

# add fossilized root states
arboreal_fossilized_relaxed <-as.factor(c(setNames("not.arboreal","ROOT"),
                                          setNames(as.character(arboreal_relaxed),names(arboreal_relaxed))))
padless_fossilized <-as.factor(c(setNames("0","ROOT"),
                                  setNames(as.character(toepads),names(toepads))))

#fit fixed-root models for Relaxed dataset

fit.fixed_dep_y_ARD_Relaxed<-fitPagel(tree_root,padless_fossilized,arboreal_fossilized_relaxed,dep.var="y",model="ARD") #relaxed
plot(fit.fixed_dep_y_ARD_Relaxed)
title("fit.fixed_dep_y_ARD_Relaxed")

fit.fixed_dep_x_ARD_Relaxed<-fitPagel(tree_root,padless_fossilized,arboreal_fossilized_relaxed,dep.var="x",model="ARD") #relaxed
plot(fit.fixed_dep_x_ARD_Relaxed)
title("fit.fixed_dep_x_ARD_Relaxed")

fit.fixed_dep_xy_ARD_Relaxed<-fitPagel(tree_root,padless_fossilized,arboreal_fossilized_relaxed,dep.var="xy",model="ARD") #relaxed
plot(fit.fixed_dep_xy_ARD_Relaxed)
title("fit.fixed_dep_xy_ARD_Relaxed")

aic_relaxed<-setNames(c(fit.fixed_dep_x_ARD_Relaxed$dependent.AIC,
                        fit.fixed_dep_y_ARD_Relaxed$dependent.AIC,
                        fit.fixed_dep_xy_ARD_Relaxed$dependent.AIC),
                      c("dependent x",
                        "dependent y","dependent x&y"))
aic_relaxed
aic.w(aic_relaxed)

##########################################################################
######                                                              ######
######   	Tests for Phylogenetic Signal with the D-statistic        ######
######                                                              ######
##########################################################################

#Toepads D-statistic
library(caper)
data_caper<-read.csv("SupplementaryTable1_EcologicalTraitDatabase_September2020.csv")
Concatenated <- data_caper$Concatenated
toepads <- data_caper$toepads
dfTOE <- data.frame(Concatenated, toepads)
phylod.toes <- phylo.d(dfTOE,tree, names.col = Concatenated, binvar = toepads)
phylod.toes
#Arboreality Binary RELAXED D-statisic
Concatenated <- data_caper$Concatenated
arboreal.relaxed.D.statistic <- data_caper$relaxed.arboreal
df.relaxed <- data.frame(Concatenated, arboreal.relaxed.D.statistic)
PhyloD.relaxed.arboreal <- phylo.d(df.relaxed, tree, names.col = Concatenated, binvar=arboreal.relaxed.D.statistic)
PhyloD.relaxed.arboreal

#Arboreality Binary STRICT D-statisic
Concatenated <- data_caper$Concatenated
arboreal.strict.D.statistic <- data_caper$strict.dataset
df.strict <- data.frame(Concatenated, arboreal.strict.D.statistic)
PhyloD.strict.arboreal <- phylo.d(df.strict, tree, names.col = Concatenated, binvar=arboreal.strict.D.statistic)
PhyloD.strict.arboreal
