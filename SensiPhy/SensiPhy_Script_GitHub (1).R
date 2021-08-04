#Investigating taxonomic and phylogenetic sensitivities using sensiPhy
library(sensiPhy)
library(ape)
library(phytools)
library(tidyverse)

## Read in Full Dataset
Strict.Full.Arbor.Data <- read.csv("FullStrict_TraitAnalysis.csv",row.names=1)
Strict.Full.Arbor <- as.factor(setNames(Strict.Full.Arbor.Data[,1],rownames(Strict.Full.Arbor.Data)))
Relaxed.Full.Arbor.Data <- read.csv("FullRelaxed_TraitAnalysis.csv",row.names=1)
Relaxed.Full.Arbor <- as.factor(setNames(Relaxed.Full.Arbor.Data[,1],rownames(Relaxed.Full.Arbor.Data)))
Full.Toepads <- as.factor(setNames(Strict.Full.Arbor.Data[,2],rownames(Strict.Full.Arbor.Data)))
FullTree <- read.nexus(file = "FullTree.tre")

## Read in Standard Dataset
Strict.Standard.Arbor.Data<-read.csv("StandardStrict_TraitAnalysis.csv",row.names=1)
Strict.Standard.Arbor<-as.factor(setNames(Strict.Standard.Arbor.Data[,1],rownames(Strict.Standard.Arbor.Data)))
Relaxed.Standard.Arbor.Data<-read.csv("StandardRelaxed_TraitAnalysis.csv",row.names=1)
Relaxed.Standard.Arbor<-as.factor(setNames(Relaxed.Standard.Arbor.Data[,1],rownames(Relaxed.Standard.Arbor.Data)))
Standard.Toepads <- as.factor(setNames(Strict.Standard.Arbor.Data[,2],rownames(Strict.Standard.Arbor.Data)))
StandardTree <- read.nexus(file = "StandardTree.tre")

## Model trait evolution to determine influential species

## Full Toepads

Influ.Full.Toepads.ARD <- influ_discrete(data = Full.Toepads, phy = FullTree,
                                         model = "ARD", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Full.Toepads.ARD.sensi <- Influ.Full.Toepads.ARD$sensi.estimates
Influ.Full.Toepads.ARD.sensi %>% filter(q12.perc > 20) ## high influence
Influ.Full.Toepads.ARD.sensi %>% filter(q21.perc > 20) ## high influence

Influ.Full.Toepads.ER <- influ_discrete(data = Full.Toepads, phy = FullTree,
                                        model = "ER", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Full.Toepads.ER.sensi <- Influ.Full.Toepads.ER$sensi.estimates
Influ.Full.Toepads.ER.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Full.Toepads.ER.sensi %>% filter(q21.perc > 20) ## no high influence

## Strict Full Arboreality

Influ.Strict.Full.Arbor.ARD <- influ_discrete(data = Strict.Full.Arbor, phy = FullTree,
                                              model = "ARD", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Strict.Full.Arbor.ARD.sensi <- Influ.Strict.Full.Arbor.ARD$sensi.estimates
Influ.Strict.Full.Arbor.ARD.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Strict.Full.Arbor.ARD.sensi %>% filter(q21.perc > 20) ## no high influence


Influ.Strict.Full.Arbor.ER <- influ_discrete(data = Strict.Full.Arbor, phy = FullTree,
                                             model = "ER", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Strict.Full.Arbor.ER.sensi <- Influ.Strict.Full.Arbor.ER$sensi.estimates
Influ.Strict.Full.Arbor.ER.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Strict.Full.Arbor.ER.sensi %>% filter(q21.perc > 20) ## no high influence

## Relaxed Full Arboreality

Influ.Relaxed.Full.Arbor.ARD <- influ_discrete(data = Relaxed.Full.Arbor, phy = FullTree,
                                               model = "ARD", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Relaxed.Full.Arbor.ARD.sensi <- Influ.Relaxed.Full.Arbor.ARD$sensi.estimates
Influ.Relaxed.Full.Arbor.ARD.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Relaxed.Full.Arbor.ARD.sensi %>% filter(q21.perc > 20) ## no high influence


Influ.Relaxed.Full.Arbor.ER <- influ_discrete(data = Relaxed.Full.Arbor, phy = FullTree,
                                              model = "ER", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Relaxed.Full.Arbor.ER.sensi <- Influ.Relaxed.Full.Arbor.ER$sensi.estimates
Influ.Relaxed.Full.Arbor.ER.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Relaxed.Full.Arbor.ER.sensi %>% filter(q21.perc > 20) ## no high influence

## Standard Toepads

Influ.Standard.Toepads.ARD <- influ_discrete(data = Standard.Toepads, phy = StandardTree,
                                             model = "ARD", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Standard.Toepads.ARD.sensi <- Influ.Standard.Toepads.ARD$sensi.estimates
Influ.Standard.Toepads.ARD.sensi %>% filter(q12.perc > 20) ## high influence
Influ.Standard.Toepads.ARD.sensi %>% filter(q21.perc > 20) ## high influence

Influ.Standard.Toepads.ER <- influ_discrete(data = Standard.Toepads, phy = StandardTree,
                                            model = "ER", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Standard.Toepads.ER.sensi <- Influ.Standard.Toepads.ER$sensi.estimates
Influ.Standard.Toepads.ER.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Standard.Toepads.ER.sensi %>% filter(q21.perc > 20) ## no high influence

## Strict Standard Arboreality

Influ.Strict.Standard.Arbor.ARD <- influ_discrete(data = Strict.Standard.Arbor, phy = StandardTree,
                                                  model = "ARD", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Strict.Standard.Arbor.ARD.sensi <- Influ.Strict.Standard.Arbor.ARD$sensi.estimates
Influ.Strict.Standard.Arbor.ARD.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Strict.Standard.Arbor.ARD.sensi %>% filter(q21.perc > 20) ## no high influence

Influ.Strict.Standard.Arbor.ER <- influ_discrete(data = Strict.Standard.Arbor, phy = StandardTree,
                                                 model = "ER", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Strict.Standard.Arbor.ER.sensi <- Influ.Strict.Standard.Arbor.ER$sensi.estimates
Influ.Strict.Standard.Arbor.ER.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Strict.Standard.Arbor.ER.sensi %>% filter(q21.perc > 20) ## no high influence

## Relaxed Standard Arboreality

Influ.Relaxed.Standard.Arbor.ARD <- influ_discrete(data = Relaxed.Standard.Arbor, phy = StandardTree,
                                                   model = "ARD", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Relaxed.Standard.Arbor.ARD.sensi <- Influ.Relaxed.Standard.Arbor.ARD$sensi.estimates ##nothing over 1.4%
Influ.Relaxed.Standard.Arbor.ARD.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Relaxed.Standard.Arbor.ARD.sensi %>% filter(q21.perc > 20) ## no high influence

Influ.Relaxed.Standard.Arbor.ER <- influ_discrete(data = Relaxed.Standard.Arbor, phy = StandardTree,
                                                  model = "ER", transform = "none", cutoff = 2, n.cores = 6, track = TRUE)

Influ.Relaxed.Standard.Arbor.ER.sensi <- Influ.Relaxed.Standard.Arbor.ER$sensi.estimates #nothing over 1%
Influ.Relaxed.Standard.Arbor.ER.sensi %>% filter(q12.perc > 20) ## no high influence
Influ.Relaxed.Standard.Arbor.ER.sensi %>% filter(q21.perc > 20) ## no high influence

## Random removal of proportion (10%, 20%, and 30%) of species

## Standard Toepads

Samp.Standard.Toepads.ER <- samp_discrete(data = Standard.Toepads,phy = StandardTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ER", n.cores = 6, track = TRUE)

Samp.Standard.Toepads.ARD <- samp_discrete(data = Standard.Toepads,phy = StandardTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ARD", n.cores = 6, track = TRUE)

## Full Toepads

Samp.Full.Toepads.ER <- samp_discrete(data = Full.Toepads,phy = FullTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ER", n.cores = 6, track = TRUE)

Samp.Full.Toepads.ARD <- samp_discrete(data = Full.Toepads,phy = FullTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ARD", n.cores = 6, track = TRUE)

## Relaxed Standard Arboreality

Samp.Relaxed.Standard.Arbor.ER <- samp_discrete(data = Relaxed.Standard.Arbor,phy = StandardTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ER", n.cores = 6, track = TRUE)

Samp.Relaxed.Standard.Arbor.ARD <- samp_discrete(data = Relaxed.Standard.Arbor,phy = StandardTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ARD", n.cores = 6, track = TRUE)

## Strict Standard Arboreality

Samp.Strict.Standard.Arbor.ER <- samp_discrete(data = Strict.Standard.Arbor,phy = StandardTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ER", n.cores = 6, track = TRUE)

Samp.Strict.Standard.Arbor.ARD <- samp_discrete(data = Strict.Standard.Arbor,phy = StandardTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ARD", n.cores = 6, track = TRUE)

## Relaxed Full Arboreality

Samp.Relaxed.Full.Arbor.ER <- samp_discrete(data = Relaxed.Full.Arbor,phy = FullTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ER", n.cores = 6, track = TRUE)

Samp.Relaxed.Full.Arbor.ARD <- samp_discrete(data = Relaxed.Full.Arbor,phy = FullTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ARD", n.cores = 6, track = TRUE)

## Strict Full Arboreality

Samp.Strict.Full.Arbor.ER <- samp_discrete(data = Strict.Full.Arbor,phy = FullTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ER", n.cores = 6, track = TRUE)

Samp.Strict.Full.Arbor.ARD <- samp_discrete(data = Strict.Full.Arbor,phy = FullTree, n.sim=25, breaks = seq(.1,.3,.1),model = "ARD", n.cores = 6, track = TRUE)

## Samp analysis

Samp.Standard.Toepads.ER.summary <- summary(Samp.Standard.Toepads.ER)
Samp.Standard.Toepads.ER.summary.q12 <- Samp.Standard.Toepads.ER.summary$`Mean q12 Change (%)`
Samp.Standard.Toepads.ER.summary.q21 <- Samp.Standard.Toepads.ER.summary$`Mean q21 Change (%)`

Samp.Standard.Toepads.ARD.summary <- summary(Samp.Standard.Toepads.ARD)
Samp.Standard.Toepads.ARD.summary.q12 <- Samp.Standard.Toepads.ARD.summary$`Mean q12 Change (%)`
Samp.Standard.Toepads.ARD.summary.q21 <- Samp.Standard.Toepads.ARD.summary$`Mean q21 Change (%)`

Samp.Full.Toepads.ER.summary <- summary(Samp.Full.Toepads.ER)
Samp.Full.Toepads.ER.summary.q12 <- Samp.Full.Toepads.ER.summary$`Mean q12 Change (%)`
Samp.Full.Toepads.ER.summary.q21 <- Samp.Full.Toepads.ER.summary$`Mean q21 Change (%)`

Samp.Full.Toepads.ARD.summary <- summary(Samp.Full.Toepads.ARD)
Samp.Full.Toepads.ARD.summary.q12 <- Samp.Full.Toepads.ARD.summary$`Mean q12 Change (%)`
Samp.Full.Toepads.ARD.summary.q21 <- Samp.Full.Toepads.ARD.summary$`Mean q21 Change (%)`

## Relaxed Standard Arboreality


Samp.Relaxed.Standard.Arbor.ER.summary <- summary(Samp.Relaxed.Standard.Arbor.ER)
Samp.Relaxed.Standard.Arbor.ER.summary.q12 <- Samp.Relaxed.Standard.Arbor.ER.summary$`Mean q12 Change (%)`
Samp.Relaxed.Standard.Arbor.ER.summary.q21 <- Samp.Relaxed.Standard.Arbor.ER.summary$`Mean q21 Change (%)`

Samp.Relaxed.Standard.Arbor.ARD.summary <- summary(Samp.Relaxed.Standard.Arbor.ARD)
Samp.Relaxed.Standard.Arbor.ARD.summary.q12 <- Samp.Relaxed.Standard.Arbor.ARD.summary$`Mean q12 Change (%)`
Samp.Relaxed.Standard.Arbor.ARD.summary.q21 <- Samp.Relaxed.Standard.Arbor.ARD.summary$`Mean q21 Change (%)`

## Strict Standard Arboreality

Samp.Strict.Standard.Arbor.ER.summary <- summary(Samp.Strict.Standard.Arbor.ER)
Samp.Strict.Standard.Arbor.ER.summary.q12 <- Samp.Strict.Standard.Arbor.ER.summary$`Mean q12 Change (%)`
Samp.Strict.Standard.Arbor.ER.summary.q21 <- Samp.Strict.Standard.Arbor.ER.summary$`Mean q21 Change (%)`

Samp.Strict.Standard.Arbor.ARD.summary <- summary(Samp.Strict.Standard.Arbor.ARD)
Samp.Strict.Standard.Arbor.ARD.summary.q12 <- Samp.Strict.Standard.Arbor.ARD.summary$`Mean q12 Change (%)`
Samp.Strict.Standard.Arbor.ARD.summary.q21 <- Samp.Strict.Standard.Arbor.ARD.summary$`Mean q21 Change (%)`

## Relaxed Full Arboreality

Samp.Relaxed.Full.Arbor.ER.summary <- summary(Samp.Relaxed.Full.Arbor.ER)
Samp.Relaxed.Full.Arbor.ER.summary.q12 <- Samp.Relaxed.Full.Arbor.ER.summary$`Mean q12 Change (%)`
Samp.Relaxed.Full.Arbor.ER.summary.q21 <- Samp.Relaxed.Full.Arbor.ER.summary$`Mean q21 Change (%)`

Samp.Relaxed.Full.Arbor.ARD.summary <- summary(Samp.Relaxed.Full.Arbor.ARD)
Samp.Relaxed.Full.Arbor.ARD.summary.q12 <- Samp.Relaxed.Full.Arbor.ARD.summary$`Mean q12 Change (%)`
Samp.Relaxed.Full.Arbor.ARD.summary.q21 <- Samp.Relaxed.Full.Arbor.ARD.summary$`Mean q21 Change (%)`

## Strict Full Arboreality

Samp.Strict.Full.Arbor.ER.summary <- summary(Samp.Strict.Full.Arbor.ER)
Samp.Strict.Full.Arbor.ER.summary.q12 <- Samp.Strict.Full.Arbor.ER.summary$`Mean q12 Change (%)`
Samp.Strict.Full.Arbor.ER.summary.q21 <- Samp.Strict.Full.Arbor.ER.summary$`Mean q21 Change (%)`

Samp.Strict.Full.Arbor.ARD.summary <- summary(Samp.Strict.Full.Arbor.ARD)
Samp.Strict.Full.Arbor.ARD.summary.q12 <- Samp.Strict.Full.Arbor.ARD.summary$`Mean q12 Change (%)`
Samp.Strict.Full.Arbor.ARD.summary.q21 <- Samp.Strict.Full.Arbor.ARD.summary$`Mean q21 Change (%)`


## Put it all together in a table and plot
Sampling.q12 <- cbind(Samp.Standard.Toepads.ER.summary.q12,
                      Samp.Standard.Toepads.ARD.summary.q12,
                      Samp.Full.Toepads.ER.summary.q12,
                      Samp.Full.Toepads.ARD.summary.q12,
                      Samp.Relaxed.Standard.Arbor.ER.summary.q12,
                      Samp.Relaxed.Standard.Arbor.ARD.summary.q12,
                      Samp.Strict.Standard.Arbor.ER.summary.q12,
                      Samp.Strict.Standard.Arbor.ARD.summary.q12,
                      Samp.Relaxed.Full.Arbor.ER.summary.q12,
                      Samp.Relaxed.Full.Arbor.ARD.summary.q12,
                      Samp.Strict.Full.Arbor.ER.summary.q12,
                      Samp.Strict.Full.Arbor.ARD.summary.q12)

Sampling.q12.Analysis <- as.data.frame(Sampling.q12, row.names = c("10% Removed",
                                                                   "20% Removed",
                                                                   "30% Removed"))

Sampling.q12.Analysis.Final <- tibble::rownames_to_column(Sampling.q12.Analysis, "LevelofChange")

#write.csv(Sampling.q12.Analysis.Final, "Sampling.q12.Analysis.Final.23Feb21.csv")

Sampling.q12.Analysis.Final.Melted <- melt((Sampling.q12.Analysis.Final))

Sampling.q21 <- cbind(Samp.Standard.Toepads.ER.summary.q21,
                      Samp.Standard.Toepads.ARD.summary.q21,
                      Samp.Full.Toepads.ER.summary.q21,
                      Samp.Full.Toepads.ARD.summary.q21,
                      Samp.Relaxed.Standard.Arbor.ER.summary.q21,
                      Samp.Relaxed.Standard.Arbor.ARD.summary.q21,
                      Samp.Strict.Standard.Arbor.ER.summary.q21,
                      Samp.Strict.Standard.Arbor.ARD.summary.q21,
                      Samp.Relaxed.Full.Arbor.ER.summary.q21,
                      Samp.Relaxed.Full.Arbor.ARD.summary.q21,
                      Samp.Strict.Full.Arbor.ER.summary.q21,
                      Samp.Strict.Full.Arbor.ARD.summary.q21)


Sampling.q21.Analysis <- as.data.frame(Sampling.q21, row.names = c("10% Removed",
                                                                   "20% Removed",
                                                                   "30% Removed"))

Sampling.q21.Analysis.Final <- tibble::rownames_to_column(Sampling.q21.Analysis, "LevelofChange")

#write.csv(Sampling.q21.Analysis.Final, "Sampling.q21.Analysis.Final.23Feb21.csv")

library(reshape2); library(Hmisc)
Sampling.q21.Analysis.Final.Melted <- melt((Sampling.q21.Analysis.Final))

Condition.Dataset.Trait.Model <- rev(c("Full Strict Arboreality ARD", 
                                       "Full Strict Arboreality ER", 
                                       "Full Relaxed Arboreality ARD",
                                       "Full Relaxed Arboreality ER",
                                       "Standard Strict Arboreality ARD",
                                       "Standard Strict Arboreality ER",
                                       "Standard Relaxed Arboreality ARD",
                                       "Standard Relaxed Arboreality ER",
                                       "Full Toepads ARD",
                                       "Full Toepads ER",
                                       "Standard Toepads ARD",
                                       "Standard Toepads ER"))

Sampling.q21.Analysis.Final.Melted.Jitter <- ggplot(Sampling.q21.Analysis.Final.Melted, aes(x=variable, y = value))+ geom_jitter(aes(color = LevelofChange), 
                                                                                                                                 position = position_jitter(0), size = 3) +
  xlab("") + labs(y=expression(Mean~"%"~Change~"in"~q["21"])) +scale_color_manual(values=c("#FFEDA0", "#FEB24C" ,"#F03B20"))+theme_classic()+coord_flip(ylim=c(0,100))+ scale_x_discrete(labels= Condition.Dataset.Trait.Model)


Sampling.q12.Analysis.Final.Melted.Jitter <- ggplot(Sampling.q12.Analysis.Final.Melted, aes(x=variable, y = value))+ geom_jitter(aes(color = LevelofChange), 
                                                                                                                                 position = position_jitter(0), size = 3) +
  xlab("") + labs(y=expression(Mean~"%"~Change~"in"~q["12"]))+scale_color_manual(values=c("#FFEDA0", "#FEB24C" ,"#F03B20"))+theme_classic()+coord_flip(ylim=c(0,100))+ scale_x_discrete(labels= Condition.Dataset.Trait.Model)



plot_grid(Sampling.q12.Analysis.Final.Melted.Jitter,Sampling.q21.Analysis.Final.Melted.Jitter, ncol = 1,align = "v",labels=c("A", "B"))

