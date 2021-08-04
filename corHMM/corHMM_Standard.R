## Standard dataset corHMM analyses

## Standard corHMM models
#Load packages
library(ape)
library(corHMM)
library(phytools)
library(geiger)

#Transition Rates Analysis corHMM
Standard_Strict_Dataset <- read.csv("StandardStrict_TraitAnalysis.csv")
Standard_Relaxed_Dataset <- read.csv("StandardRelaxed_TraitAnalysis.csv")
StandardTree <- read.nexus(file = "StandardTree.tre")

##STRICT with no hiddden classes###

#ARD No Dual Transitions
Strict_ARD_NoDual <- getStateMat4Dat(Standard_Strict_Dataset, model="ARD",dual = FALSE)
Strict_ARD_NoDual_NoHidden_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ARD_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ARD Dual Transitions Allowed
Strict_ARD_Dual <- getStateMat4Dat(Standard_Strict_Dataset, model="ARD", dual = TRUE)
Strict_ARD_Dual_NoHidden_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ARD_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER No Dual Transitions
Strict_ER_NoDual <- getStateMat4Dat(Standard_Strict_Dataset, model="ER", dual = FALSE)
Strict_ER_NoDual_NoHidden_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ER_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER Dual Transitions Allowed
Strict_ER_Dual <- getStateMat4Dat(Standard_Strict_Dataset, model="ER", dual = TRUE)
Strict_ER_Dual_NoHidden_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ER_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with no hiddden classes###

#ARD No Dual Transitions
Relaxed_ARD_NoDual <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ARD", dual = FALSE)
Relaxed_ARD_NoDual_NoHidden_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ARD_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ARD Dual Transitions Allowed
Relaxed_ARD_Dual <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ARD", dual = TRUE)
Relaxed_ARD_Dual_NoHidden_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ARD_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER No Dual Transitions
Relaxed_ER_NoDual <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ER", dual = FALSE)
Relaxed_ER_NoDual_NoHidden_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ER_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER Dual Transitions Allowed
Relaxed_ER_Dual <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ER", dual = TRUE)
Relaxed_ER_Dual_NoHidden_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ER_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###STRICT with rate.cat = 2 ###

#ARD No Dual Transitions
ARD_strict_NoDual_Hidden2_R1 <- getStateMat4Dat(Standard_Strict_Dataset, model="ARD",dual = FALSE)$rate.mat
ARD_strict_NoDual_Hidden2_R2 <- getStateMat4Dat(Standard_Strict_Dataset, model="ARD",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ARD_NoDual_HiddenRatCat2 <- getFullMat(list(ARD_strict_NoDual_Hidden2_R1, ARD_strict_NoDual_Hidden2_R2), Pp2)
Strict_ARD_NoDual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ARD_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ARD Dual Transitions Allowed
ARD_strict_Dual_Hidden2_R1 <- getStateMat4Dat(Standard_Strict_Dataset, model="ARD",dual = TRUE)$rate.mat
ARD_strict_Dual_Hidden2_R2 <- getStateMat4Dat(Standard_Strict_Dataset, model="ARD",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ARD_DualAllowed_HiddenRatCat2 <- getFullMat(list(ARD_strict_Dual_Hidden2_R1, ARD_strict_Dual_Hidden2_R2), Pp2)
Strict_ARD_DualAllowed_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ARD_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ER No Dual Transitions
ER_strict_NoDual_Hidden2_R1 <- getStateMat4Dat(Standard_Strict_Dataset, model="ER",dual = FALSE)$rate.mat
ER_strict_NoDual_Hidden2_R2 <- getStateMat4Dat(Standard_Strict_Dataset, model="ER",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ER_NoDual_HiddenRatCat2 <- getFullMat(list(ER_strict_NoDual_Hidden2_R1, ER_strict_NoDual_Hidden2_R2), Pp2)
Strict_ER_NoDual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ER_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ER Dual Transitions Allowed
ER_strict_Dual_Hidden2_R1 <- getStateMat4Dat(Standard_Strict_Dataset, model="ER",dual = TRUE)$rate.mat
ER_strict_Dual_Hidden2_R2 <- getStateMat4Dat(Standard_Strict_Dataset, model="ER",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ER_DualAllowed_HiddenRatCat2 <- getFullMat(list(ER_strict_Dual_Hidden2_R1, ER_strict_Dual_Hidden2_R2), Pp2)
Strict_ER_DualAllowed_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_ER_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with rate.cat = 2###

#ARD No Dual Transitions
ARD_relaxed_NoDual_Hidden2_R1 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ARD",dual = FALSE)$rate.mat
ARD_relaxed_NoDual_Hidden2_R2 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ARD",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ARD_NoDual_HiddenRatCat2 <- getFullMat(list(ARD_relaxed_NoDual_Hidden2_R1, ARD_relaxed_NoDual_Hidden2_R2), Pp2)
Relaxed_ARD_NoDual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ARD_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ARD Dual Transitions Allowed
ARD_relaxed_Dual_Hidden2_R1 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ARD",dual = TRUE)$rate.mat
ARD_relaxed_Dual_Hidden2_R2 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ARD",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ARD_Dual_HiddenRatCat2 <- getFullMat(list(ARD_relaxed_Dual_Hidden2_R1, ARD_relaxed_Dual_Hidden2_R2), Pp2)
Relaxed_ARD_Dual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ARD_Dual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
Relaxed_ARD_DualAllowed_HiddenRatCat2_Model <- Relaxed_ARD_Dual_HiddenRatCat2_Model
#ER No Dual Transitions
ER_relaxed_NoDual_Hidden2_R1 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ER",dual = FALSE)$rate.mat
ER_relaxed_NoDual_Hidden2_R2 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ER",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ER_NoDual_HiddenRatCat2 <- getFullMat(list(ER_relaxed_NoDual_Hidden2_R1, ER_relaxed_NoDual_Hidden2_R2), Pp2)
Relaxed_ER_NoDual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ER_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ER Dual Transitions Allowed
ER_relaxed_Dual_Hidden2_R1 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ER",dual = TRUE)$rate.mat
ER_relaxed_Dual_Hidden2_R2 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="ER",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ER_DualAllowed_HiddenRatCat2 <- getFullMat(list(ER_relaxed_Dual_Hidden2_R1, ER_relaxed_Dual_Hidden2_R2), Pp2)
Relaxed_ER_DualAllowed_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_ER_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)


##SYMMETRIC RATES MODELS

#SYM No Dual Transitions
Strict_SYM_NoDual <- getStateMat4Dat(Standard_Strict_Dataset, model="SYM",dual = FALSE)
Strict_SYM_NoDual_NoHidden_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_SYM_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#SYM Dual Transitions Allowed
Strict_SYM_Dual <- getStateMat4Dat(Standard_Strict_Dataset, model="SYM", dual = TRUE)
Strict_SYM_Dual_NoHidden_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_SYM_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with no hiddden classes###

#SYM No Dual Transitions
Relaxed_SYM_NoDual <- getStateMat4Dat(Standard_Relaxed_Dataset, model="SYM", dual = FALSE)
Relaxed_SYM_NoDual_NoHidden_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_SYM_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#SYM Dual Transitions Allowed
Relaxed_SYM_Dual <- getStateMat4Dat(Standard_Relaxed_Dataset, model="SYM", dual = TRUE)
Relaxed_SYM_Dual_NoHidden_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_SYM_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###STRICT with rate.cat = 2 ###

#SYM No Dual Transitions
SYM_strict_NoDual_Hidden2_R1 <- getStateMat4Dat(Standard_Strict_Dataset, model="SYM",dual = FALSE)$rate.mat
SYM_strict_NoDual_Hidden2_R2 <- getStateMat4Dat(Standard_Strict_Dataset, model="SYM",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_SYM_NoDual_HiddenRatCat2 <- getFullMat(list(SYM_strict_NoDual_Hidden2_R1, SYM_strict_NoDual_Hidden2_R2), Pp2)
Strict_SYM_NoDual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_SYM_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#SYM Dual Transitions Allowed
SYM_strict_Dual_Hidden2_R1 <- getStateMat4Dat(Standard_Strict_Dataset, model="SYM",dual = TRUE)$rate.mat
SYM_strict_Dual_Hidden2_R2 <- getStateMat4Dat(Standard_Strict_Dataset, model="SYM",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_SYM_DualAllowed_HiddenRatCat2 <- getFullMat(list(SYM_strict_Dual_Hidden2_R1, SYM_strict_Dual_Hidden2_R2), Pp2)
Strict_SYM_DualAllowed_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Strict_Dataset,rate.mat = Strict_SYM_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with rate.cat = 2###

#SYM No Dual Transitions
SYM_relaxed_NoDual_Hidden2_R1 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="SYM",dual = FALSE)$rate.mat
SYM_relaxed_NoDual_Hidden2_R2 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="SYM",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_SYM_NoDual_HiddenRatCat2 <- getFullMat(list(SYM_relaxed_NoDual_Hidden2_R1, SYM_relaxed_NoDual_Hidden2_R2), Pp2)
Relaxed_SYM_NoDual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_SYM_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#SYM Dual Transitions Allowed
SYM_relaxed_Dual_Hidden2_R1 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="SYM",dual = TRUE)$rate.mat
SYM_relaxed_Dual_Hidden2_R2 <- getStateMat4Dat(Standard_Relaxed_Dataset, model="SYM",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_SYM_Dual_HiddenRatCat2 <- getFullMat(list(SYM_relaxed_Dual_Hidden2_R1, SYM_relaxed_Dual_Hidden2_R2), Pp2)
Relaxed_SYM_Dual_HiddenRatCat2_Model <- corHMM(StandardTree,Standard_Relaxed_Dataset,rate.mat = Relaxed_SYM_Dual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
Relaxed_SYM_DualAllowed_HiddenRatCat2_Model <- Relaxed_SYM_Dual_HiddenRatCat2_Model

## Standard Strict table

AICc_Strict_Standard <- c(Strict_ER_Dual_NoHidden_Model$AICc, 
                          Strict_ER_NoDual_NoHidden_Model$AICc,
                          Strict_ER_DualAllowed_HiddenRatCat2_Model$AICc, 
                          Strict_ER_NoDual_HiddenRatCat2_Model$AICc, 
                          Strict_ARD_DualAllowed_HiddenRatCat2_Model$AICc, 
                          Strict_ARD_Dual_NoHidden_Model$AICc, 
                          Strict_ARD_NoDual_HiddenRatCat2_Model$AICc, 
                          Strict_ARD_NoDual_NoHidden_Model$AICc,
                          Strict_SYM_DualAllowed_HiddenRatCat2_Model$AICc,
                          Strict_SYM_Dual_NoHidden_Model$AICc,
                          Strict_SYM_NoDual_HiddenRatCat2_Model$AICc,
                          Strict_SYM_NoDual_NoHidden_Model$AICc)

LogL_Strict_Standard <- c(Strict_ER_Dual_NoHidden_Model$loglik, 
                          Strict_ER_NoDual_NoHidden_Model$loglik,
                          Strict_ER_DualAllowed_HiddenRatCat2_Model$loglik, 
                          Strict_ER_NoDual_HiddenRatCat2_Model$loglik, 
                          Strict_ARD_DualAllowed_HiddenRatCat2_Model$loglik, 
                          Strict_ARD_Dual_NoHidden_Model$loglik, 
                          Strict_ARD_NoDual_HiddenRatCat2_Model$loglik, 
                          Strict_ARD_NoDual_NoHidden_Model$loglik,
                          Strict_SYM_DualAllowed_HiddenRatCat2_Model$loglik,
                          Strict_SYM_Dual_NoHidden_Model$loglik,
                          Strict_SYM_NoDual_HiddenRatCat2_Model$loglik,
                          Strict_SYM_NoDual_NoHidden_Model$loglik)


All_Strict_AICc_Standard <-as.data.frame(AICc_Strict_Standard, row.names = c("Strict_ER_Dual_NoHidden_Model", 
                                                                             "Strict_ER_NoDual_NoHidden_Model",
                                                                             "Strict_ER_DualAllowed_HiddenRatCat2_Model", 
                                                                             "Strict_ER_NoDual_HiddenRatCat2_Model", 
                                                                             "Strict_ARD_DualAllowed_HiddenRatCat2_Model", 
                                                                             "Strict_ARD_Dual_NoHidden_Model", 
                                                                             "Strict_ARD_NoDual_HiddenRatCat2_Model", 
                                                                             "Strict_ARD_NoDual_NoHidden_Model",
                                                                             "Strict_SYM_DualAllowed_HiddenRatCat2_Model",
                                                                             "Strict_SYM_Dual_NoHidden_Model",
                                                                             "Strict_SYM_NoDual_HiddenRatCat2_Model",
                                                                             "Strict_SYM_NoDual_NoHidden_Model"))

All_Strict_Standard <-as.data.frame(LogL_Strict_Standard, row.names = c("Strict_ER_Dual_NoHidden_Model", 
                                                                        "Strict_ER_NoDual_NoHidden_Model",
                                                                        "Strict_ER_DualAllowed_HiddenRatCat2_Model", 
                                                                        "Strict_ER_NoDual_HiddenRatCat2_Model", 
                                                                        "Strict_ARD_DualAllowed_HiddenRatCat2_Model", 
                                                                        "Strict_ARD_Dual_NoHidden_Model", 
                                                                        "Strict_ARD_NoDual_HiddenRatCat2_Model", 
                                                                        "Strict_ARD_NoDual_NoHidden_Model",
                                                                        "Strict_SYM_DualAllowed_HiddenRatCat2_Model",
                                                                        "Strict_SYM_Dual_NoHidden_Model",
                                                                        "Strict_SYM_NoDual_HiddenRatCat2_Model",
                                                                        "Strict_SYM_NoDual_NoHidden_Model"))


All_Strict_Standard$AICc<-All_Strict_AICc_Standard #Strict_SYM_DualAllowed_HiddenRatCat2_Model wins at 98.9% wAIC
aicw_strict_Standard <- aicw(AICc_Strict_Standard)
All_Strict_Standard$AICw<-aicw_strict_Standard$w


Models_StandardStrict <- c(list(Strict_ARD_Dual_NoHidden_Model,
                 Strict_ARD_DualAllowed_HiddenRatCat2_Model,
                 Strict_ARD_NoDual_HiddenRatCat2_Model,
                 Strict_ARD_NoDual_NoHidden_Model,
                 Strict_ER_Dual_NoHidden_Model,
                 Strict_ER_DualAllowed_HiddenRatCat2_Model,
                 Strict_ER_NoDual_HiddenRatCat2_Model,
                 Strict_ER_NoDual_NoHidden_Model,
                 Strict_SYM_Dual_NoHidden_Model,
                 Strict_SYM_DualAllowed_HiddenRatCat2_Model,
                 Strict_SYM_NoDual_HiddenRatCat2_Model,
                 Strict_SYM_NoDual_NoHidden_Model))

AllTable_StandardStrict <- matrix(0, length(Models_StandardStrict), 4)
rownames(AllTable_StandardStrict) <- rep(c("Standard_Strict_ARD_Dual_NoHidden_Model",
                            "Standard_Strict_ARD_DualAllowed_HiddenRatCat2_Model",
                            "Standard_Strict_ARD_NoDual_HiddenRatCat2_Model",
                            "Standard_Strict_ARD_NoDual_NoHidden_Model",
                            "Standard_Strict_ER_Dual_NoHidden_Model",
                            "Standard_Strict_ER_DualAllowed_HiddenRatCat2_Model",
                            "Standard_Strict_ER_NoDual_HiddenRatCat2_Model",
                            "Standard_Strict_ER_NoDual_NoHidden_Model",
                            "Standard_Strict_SYM_Dual_NoHidden_Model",
                            "Standard_Strict_SYM_DualAllowed_HiddenRatCat2_Model",
                            "Standard_Strict_SYM_NoDual_HiddenRatCat2_Model",
                            "Standard_Strict_SYM_NoDual_NoHidden_Model"))

colnames(AllTable_StandardStrict) <- c("k.rate", "AICc", "MeanRate", "LogL")

AllTable_StandardStrict <- as.data.frame(AllTable_StandardStrict)

#For loop to put everything into a nice table
for(i in 1:length(Models_StandardStrict)){
  k.rate <- max(Models_StandardStrict[[i]]$index.mat, na.rm = TRUE)
  AICc <- round(Models_StandardStrict[[i]]$AICc, 2)
  MeanRate <- round(mean(Models_StandardStrict[[i]]$solution, na.rm = TRUE),2)
  AllTable_StandardStrict$k.rate[i] <- print(k.rate)
  AllTable_StandardStrict$AICc[i] <- print(AICc)
  AllTable_StandardStrict_AICw <- aicw(AllTable_StandardStrict$AICc)
  AllTable_StandardStrict$AICw <- round(AllTable_StandardStrict_AICw$w,2)
  AllTable_StandardStrict$LogL <- All_Strict_Standard$LogL_Strict_Standard
  AllTable_StandardStrict$MeanRate[i] <- print(MeanRate)
  View(AllTable_StandardStrict)
}


## Standard Relaxed table

AICc_Relaxed_Standard <- c(Relaxed_ER_Dual_NoHidden_Model$AICc, 
                           Relaxed_ER_NoDual_NoHidden_Model$AICc,
                           Relaxed_ER_DualAllowed_HiddenRatCat2_Model$AICc, 
                           Relaxed_ER_NoDual_HiddenRatCat2_Model$AICc, 
                           Relaxed_ARD_Dual_HiddenRatCat2_Model$AICc, 
                           Relaxed_ARD_Dual_NoHidden_Model$AICc, 
                           Relaxed_ARD_NoDual_HiddenRatCat2_Model$AICc, 
                           Relaxed_ARD_NoDual_NoHidden_Model$AICc,
                           Relaxed_SYM_Dual_HiddenRatCat2_Model$AICc,
                           Relaxed_SYM_Dual_NoHidden_Model$AICc,
                           Relaxed_SYM_NoDual_HiddenRatCat2_Model$AICc,
                           Relaxed_SYM_NoDual_NoHidden_Model$AICc)

LogL_Relaxed_Standard <- c(Relaxed_ER_Dual_NoHidden_Model$loglik, 
                           Relaxed_ER_NoDual_NoHidden_Model$loglik,
                           Relaxed_ER_DualAllowed_HiddenRatCat2_Model$loglik, 
                           Relaxed_ER_NoDual_HiddenRatCat2_Model$loglik, 
                           Relaxed_ARD_Dual_HiddenRatCat2_Model$loglik, 
                           Relaxed_ARD_Dual_NoHidden_Model$loglik, 
                           Relaxed_ARD_NoDual_HiddenRatCat2_Model$loglik, 
                           Relaxed_ARD_NoDual_NoHidden_Model$loglik,
                           Relaxed_SYM_Dual_HiddenRatCat2_Model$loglik,
                           Relaxed_SYM_Dual_NoHidden_Model$loglik,
                           Relaxed_SYM_NoDual_HiddenRatCat2_Model$loglik,
                           Relaxed_SYM_NoDual_NoHidden_Model$loglik)


All_Relaxed_AICc_Standard <-as.data.frame(AICc_Relaxed_Standard, row.names = c("Relaxed_ER_Dual_NoHidden_Model", 
                                                                               "Relaxed_ER_NoDual_NoHidden_Model",
                                                                               "Relaxed_ER_DualAllowed_HiddenRatCat2_Model", 
                                                                               "Relaxed_ER_NoDual_HiddenRatCat2_Model", 
                                                                               "Relaxed_ARD_DualAllowed_HiddenRatCat2_Model", 
                                                                               "Relaxed_ARD_Dual_NoHidden_Model", 
                                                                               "Relaxed_ARD_NoDual_HiddenRatCat2_Model", 
                                                                               "Relaxed_ARD_NoDual_NoHidden_Model",
                                                                               "Relaxed_SYM_DualAllowed_HiddenRatCat2_Model",
                                                                               "Relaxed_SYM_Dual_NoHidden_Model",
                                                                               "Relaxed_SYM_NoDual_HiddenRatCat2_Model",
                                                                               "Relaxed_SYM_NoDual_NoHidden_Model"))

All_Relaxed_Standard <-as.data.frame(LogL_Relaxed_Standard, row.names = c("Relaxed_ER_Dual_NoHidden_Model", 
                                                                          "Relaxed_ER_NoDual_NoHidden_Model",
                                                                          "Relaxed_ER_DualAllowed_HiddenRatCat2_Model", 
                                                                          "Relaxed_ER_NoDual_HiddenRatCat2_Model", 
                                                                          "Relaxed_ARD_DualAllowed_HiddenRatCat2_Model", 
                                                                          "Relaxed_ARD_Dual_NoHidden_Model", 
                                                                          "Relaxed_ARD_NoDual_HiddenRatCat2_Model", 
                                                                          "Relaxed_ARD_NoDual_NoHidden_Model",
                                                                          "Relaxed_SYM_DualAllowed_HiddenRatCat2_Model",
                                                                          "Relaxed_SYM_Dual_NoHidden_Model",
                                                                          "Relaxed_SYM_NoDual_HiddenRatCat2_Model",
                                                                          "Relaxed_SYM_NoDual_NoHidden_Model"))


All_Relaxed_Standard$AICc<-All_Relaxed_AICc_Standard #Relaxed_SYM_DualAllowed_HiddenRatCat2_Model wins at 98.9% wAIC
aicw_Relaxed_Standard <- aicw(AICc_Relaxed_Standard)
All_Relaxed_Standard$AICw<-aicw_Relaxed_Standard$w


Models_StandardRelaxed <- c(list(Relaxed_ARD_Dual_NoHidden_Model,
                                 Relaxed_ARD_Dual_HiddenRatCat2_Model,
                                 Relaxed_ARD_NoDual_HiddenRatCat2_Model,
                                 Relaxed_ARD_NoDual_NoHidden_Model,
                                 Relaxed_ER_Dual_NoHidden_Model,
                                 Relaxed_ER_DualAllowed_HiddenRatCat2_Model,
                                 Relaxed_ER_NoDual_HiddenRatCat2_Model,
                                 Relaxed_ER_NoDual_NoHidden_Model,
                                 Relaxed_SYM_Dual_NoHidden_Model,
                                 Relaxed_SYM_Dual_HiddenRatCat2_Model,
                                 Relaxed_SYM_NoDual_HiddenRatCat2_Model,
                                 Relaxed_SYM_NoDual_NoHidden_Model))

AllTable_StandardRelaxed <- matrix(0, length(Models_StandardRelaxed), 4)
rownames(AllTable_StandardRelaxed) <- rep(c("Standard_Relaxed_ARD_Dual_NoHidden_Model",
                                            "Standard_Relaxed_ARD_DualAllowed_HiddenRatCat2_Model",
                                            "Standard_Relaxed_ARD_NoDual_HiddenRatCat2_Model",
                                            "Standard_Relaxed_ARD_NoDual_NoHidden_Model",
                                            "Standard_Relaxed_ER_Dual_NoHidden_Model",
                                            "Standard_Relaxed_ER_DualAllowed_HiddenRatCat2_Model",
                                            "Standard_Relaxed_ER_NoDual_HiddenRatCat2_Model",
                                            "Standard_Relaxed_ER_NoDual_NoHidden_Model",
                                            "Standard_Relaxed_SYM_Dual_NoHidden_Model",
                                            "Standard_Relaxed_SYM_DualAllowed_HiddenRatCat2_Model",
                                            "Standard_Relaxed_SYM_NoDual_HiddenRatCat2_Model",
                                            "Standard_Relaxed_SYM_NoDual_NoHidden_Model"))

colnames(AllTable_StandardRelaxed) <- c("k.rate", "AICc", "MeanRate", "LogL")

AllTable_StandardRelaxed <- as.data.frame(AllTable_StandardRelaxed)

#For loop to put everything together in a table
for(i in 1:length(Models_StandardRelaxed)){
  k.rate <- max(Models_StandardRelaxed[[i]]$index.mat, na.rm = TRUE)
  AICc <- round(Models_StandardRelaxed[[i]]$AICc, 2)
  MeanRate <- round(mean(Models_StandardRelaxed[[i]]$solution, na.rm = TRUE),2)
  AllTable_StandardRelaxed$k.rate[i] <- print(k.rate)
  AllTable_StandardRelaxed$AICc[i] <- print(AICc)
  AllTable_StandardRelaxed_Print_AICw <- aicw(AllTable_StandardRelaxed$AICc)
  AllTable_StandardRelaxed$AICw <- round(AllTable_StandardRelaxed_Print_AICw$w,2)
  AllTable_StandardRelaxed$LogL <- All_Relaxed_Standard$LogL_Relaxed_Standard
  AllTable_StandardRelaxed$MeanRate[i] <- print(MeanRate)
  View(AllTable_StandardRelaxed)
}


