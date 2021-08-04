## Full dataset corHMM analyses

## Full corHMM models
#Load packages
library(ape)
library(corHMM)
library(phytools)
library(geiger)

#Transition Rates Analysis corHMM
Full_Strict_Dataset <- read.csv("FullStrict_TraitAnalysis.csv")
Full_Relaxed_Dataset <- read.csv("FullRelaxed_TraitAnalysis.csv")
tree <- read.nexus(file = "FullTree.tre")

##STRICT with no hiddden classes###

#ARD No Dual Transitions
Strict_ARD_NoDual <- getStateMat4Dat(Full_Strict_Dataset, model="ARD",dual = FALSE)
Strict_ARD_NoDual_NoHidden_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ARD_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ARD Dual Transitions Allowed
Strict_ARD_Dual <- getStateMat4Dat(Full_Strict_Dataset, model="ARD", dual = TRUE)
Strict_ARD_Dual_NoHidden_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ARD_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER No Dual Transitions
Strict_ER_NoDual <- getStateMat4Dat(Full_Strict_Dataset, model="ER", dual = FALSE)
Strict_ER_NoDual_NoHidden_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ER_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER Dual Transitions Allowed
Strict_ER_Dual <- getStateMat4Dat(Full_Strict_Dataset, model="ER", dual = TRUE)
Strict_ER_Dual_NoHidden_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ER_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with no hiddden classes###

#ARD No Dual Transitions
Relaxed_ARD_NoDual <- getStateMat4Dat(Full_Relaxed_Dataset, model="ARD", dual = FALSE)
Relaxed_ARD_NoDual_NoHidden_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ARD_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ARD Dual Transitions Allowed
Relaxed_ARD_Dual <- getStateMat4Dat(Full_Relaxed_Dataset, model="ARD", dual = TRUE)
Relaxed_ARD_Dual_NoHidden_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ARD_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER No Dual Transitions
Relaxed_ER_NoDual <- getStateMat4Dat(Full_Relaxed_Dataset, model="ER", dual = FALSE)
Relaxed_ER_NoDual_NoHidden_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ER_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#ER Dual Transitions Allowed
Relaxed_ER_Dual <- getStateMat4Dat(Full_Relaxed_Dataset, model="ER", dual = TRUE)
Relaxed_ER_Dual_NoHidden_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ER_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###STRICT with rate.cat = 2 ###

#ARD No Dual Transitions
ARD_strict_NoDual_Hidden2_R1 <- getStateMat4Dat(Full_Strict_Dataset, model="ARD",dual = FALSE)$rate.mat
ARD_strict_NoDual_Hidden2_R2 <- getStateMat4Dat(Full_Strict_Dataset, model="ARD",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ARD_NoDual_HiddenRatCat2 <- getFullMat(list(ARD_strict_NoDual_Hidden2_R1, ARD_strict_NoDual_Hidden2_R2), Pp2)
Strict_ARD_NoDual_HiddenRatCat2_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ARD_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ARD Dual Transitions Allowed
ARD_strict_Dual_Hidden2_R1 <- getStateMat4Dat(Full_Strict_Dataset, model="ARD",dual = TRUE)$rate.mat
ARD_strict_Dual_Hidden2_R2 <- getStateMat4Dat(Full_Strict_Dataset, model="ARD",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ARD_DualAllowed_HiddenRatCat2 <- getFullMat(list(ARD_strict_Dual_Hidden2_R1, ARD_strict_Dual_Hidden2_R2), Pp2)
Strict_ARD_DualAllowed_HiddenRatCat2_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ARD_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ER No Dual Transitions
ER_strict_NoDual_Hidden2_R1 <- getStateMat4Dat(Full_Strict_Dataset, model="ER",dual = FALSE)$rate.mat
ER_strict_NoDual_Hidden2_R2 <- getStateMat4Dat(Full_Strict_Dataset, model="ER",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ER_NoDual_HiddenRatCat2 <- getFullMat(list(ER_strict_NoDual_Hidden2_R1, ER_strict_NoDual_Hidden2_R2), Pp2)
Strict_ER_NoDual_HiddenRatCat2_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ER_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ER Dual Transitions Allowed
ER_strict_Dual_Hidden2_R1 <- getStateMat4Dat(Full_Strict_Dataset, model="ER",dual = TRUE)$rate.mat
ER_strict_Dual_Hidden2_R2 <- getStateMat4Dat(Full_Strict_Dataset, model="ER",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_ER_DualAllowed_HiddenRatCat2 <- getFullMat(list(ER_strict_Dual_Hidden2_R1, ER_strict_Dual_Hidden2_R2), Pp2)
Strict_ER_DualAllowed_HiddenRatCat2_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_ER_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with rate.cat = 2###

#ARD No Dual Transitions
ARD_relaxed_NoDual_Hidden2_R1 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ARD",dual = FALSE)$rate.mat
ARD_relaxed_NoDual_Hidden2_R2 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ARD",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ARD_NoDual_HiddenRatCat2 <- getFullMat(list(ARD_relaxed_NoDual_Hidden2_R1, ARD_relaxed_NoDual_Hidden2_R2), Pp2)
Relaxed_ARD_NoDual_HiddenRatCat2_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ARD_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ARD Dual Transitions Allowed
ARD_relaxed_Dual_Hidden2_R1 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ARD",dual = TRUE)$rate.mat
ARD_relaxed_Dual_Hidden2_R2 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ARD",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ARD_Dual_HiddenRatCat2 <- getFullMat(list(ARD_relaxed_Dual_Hidden2_R1, ARD_relaxed_Dual_Hidden2_R2), Pp2)
Relaxed_ARD_Dual_HiddenRatCat2_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ARD_Dual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ER No Dual Transitions
ER_relaxed_NoDual_Hidden2_R1 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ER",dual = FALSE)$rate.mat
ER_relaxed_NoDual_Hidden2_R2 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ER",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ER_NoDual_HiddenRatCat2 <- getFullMat(list(ER_relaxed_NoDual_Hidden2_R1, ER_relaxed_NoDual_Hidden2_R2), Pp2)
Relaxed_ER_NoDual_HiddenRatCat2_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ER_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#ER Dual Transitions Allowed
ER_relaxed_Dual_Hidden2_R1 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ER",dual = TRUE)$rate.mat
ER_relaxed_Dual_Hidden2_R2 <- getStateMat4Dat(Full_Relaxed_Dataset, model="ER",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_ER_DualAllowed_HiddenRatCat2 <- getFullMat(list(ER_relaxed_Dual_Hidden2_R1, ER_relaxed_Dual_Hidden2_R2), Pp2)
Relaxed_ER_DualAllowed_HiddenRatCat2_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_ER_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)


##SYMMETRIC RATES MODELS

#SYM No Dual Transitions
Strict_SYM_NoDual <- getStateMat4Dat(Full_Strict_Dataset, model="SYM",dual = FALSE)
Strict_SYM_NoDual_NoHidden_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_SYM_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#SYM Dual Transitions Allowed
Strict_SYM_Dual <- getStateMat4Dat(Full_Strict_Dataset, model="SYM", dual = TRUE)
Strict_SYM_Dual_NoHidden_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_SYM_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with no hiddden classes###

#SYM No Dual Transitions
Relaxed_SYM_NoDual <- getStateMat4Dat(Full_Relaxed_Dataset, model="SYM", dual = FALSE)
Relaxed_SYM_NoDual_NoHidden_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_SYM_NoDual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)
#SYM Dual Transitions Allowed
Relaxed_SYM_Dual <- getStateMat4Dat(Full_Relaxed_Dataset, model="SYM", dual = TRUE)
Relaxed_SYM_Dual_NoHidden_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_SYM_Dual$rate.mat, rate.cat=1,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###STRICT with rate.cat = 2 ###

#SYM No Dual Transitions
SYM_strict_NoDual_Hidden2_R1 <- getStateMat4Dat(Full_Strict_Dataset, model="SYM",dual = FALSE)$rate.mat
SYM_strict_NoDual_Hidden2_R2 <- getStateMat4Dat(Full_Strict_Dataset, model="SYM",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_SYM_NoDual_HiddenRatCat2 <- getFullMat(list(SYM_strict_NoDual_Hidden2_R1, SYM_strict_NoDual_Hidden2_R2), Pp2)
Strict_SYM_NoDual_HiddenRatCat2_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_SYM_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#SYM Dual Transitions Allowed
SYM_strict_Dual_Hidden2_R1 <- getStateMat4Dat(Full_Strict_Dataset, model="SYM",dual = TRUE)$rate.mat
SYM_strict_Dual_Hidden2_R2 <- getStateMat4Dat(Full_Strict_Dataset, model="SYM",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Strict_SYM_DualAllowed_HiddenRatCat2 <- getFullMat(list(SYM_strict_Dual_Hidden2_R1, SYM_strict_Dual_Hidden2_R2), Pp2)
Strict_SYM_DualAllowed_HiddenRatCat2_Model <- corHMM(tree,Full_Strict_Dataset,rate.mat = Strict_SYM_DualAllowed_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

###RELAXED with rate.cat = 2###

#SYM No Dual Transitions
SYM_relaxed_NoDual_Hidden2_R1 <- getStateMat4Dat(Full_Relaxed_Dataset, model="SYM",dual = FALSE)$rate.mat
SYM_relaxed_NoDual_Hidden2_R2 <- getStateMat4Dat(Full_Relaxed_Dataset, model="SYM",dual = FALSE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_SYM_NoDual_HiddenRatCat2 <- getFullMat(list(SYM_relaxed_NoDual_Hidden2_R1, SYM_relaxed_NoDual_Hidden2_R2), Pp2)
Relaxed_SYM_NoDual_HiddenRatCat2_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_SYM_NoDual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

#SYM Dual Transitions Allowed
SYM_relaxed_Dual_Hidden2_R1 <- getStateMat4Dat(Full_Relaxed_Dataset, model="SYM",dual = TRUE)$rate.mat
SYM_relaxed_Dual_Hidden2_R2 <- getStateMat4Dat(Full_Relaxed_Dataset, model="SYM",dual = TRUE)$rate.mat
Pp2 <- getRateCatMat(2)
Relaxed_SYM_Dual_HiddenRatCat2 <- getFullMat(list(SYM_relaxed_Dual_Hidden2_R1, SYM_relaxed_Dual_Hidden2_R2), Pp2)
Relaxed_SYM_Dual_HiddenRatCat2_Model <- corHMM(tree,Full_Relaxed_Dataset,rate.mat = Relaxed_SYM_Dual_HiddenRatCat2, rate.cat=2,root.p=c(1,0,0,0),nstarts=100,node.states="marginal",ip=1,get.tip.states = TRUE)

##Full Relaxed table

AICc_Relaxed_Full <- c(Relaxed_ER_Dual_NoHidden_Model$AICc, 
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

LogL_Relaxed_Full <- c(Relaxed_ER_Dual_NoHidden_Model$loglik, 
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


All_Relaxed_AICc_Full <-as.data.frame(AICc_Relaxed_Full, row.names = c("Relaxed_ER_Dual_NoHidden_Model", 
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

All_Relaxed_Full <-as.data.frame(LogL_Relaxed_Full, row.names = c("Relaxed_ER_Dual_NoHidden_Model", 
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


All_Relaxed_Full$AICc<-All_Relaxed_AICc_Full #Relaxed_SYM_DualAllowed_HiddenRatCat2_Model wins at 98.9% wAIC


Models_FullRelaxed <- c(list(Relaxed_ARD_Dual_NoHidden_Model,
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

AllTable_FullRelaxed <- matrix(0, length(Models_FullRelaxed), 4)
rownames(AllTable_FullRelaxed) <- rep(c("Full_Relaxed_ARD_Dual_NoHidden_Model",
                                        "Full_Relaxed_ARD_DualAllowed_HiddenRatCat2_Model",
                                        "Full_Relaxed_ARD_NoDual_HiddenRatCat2_Model",
                                        "Full_Relaxed_ARD_NoDual_NoHidden_Model",
                                        "Full_Relaxed_ER_Dual_NoHidden_Model",
                                        "Full_Relaxed_ER_DualAllowed_HiddenRatCat2_Model",
                                        "Full_Relaxed_ER_NoDual_HiddenRatCat2_Model",
                                        "Full_Relaxed_ER_NoDual_NoHidden_Model",
                                        "Full_Relaxed_SYM_Dual_NoHidden_Model",
                                        "Full_Relaxed_SYM_DualAllowed_HiddenRatCat2_Model",
                                        "Full_Relaxed_SYM_NoDual_HiddenRatCat2_Model",
                                        "Full_Relaxed_SYM_NoDual_NoHidden_Model"))

colnames(AllTable_FullRelaxed) <- c("k.rate", "AICc", "MeanRate", "LogL")

AllTable_FullRelaxed <- as.data.frame(AllTable_FullRelaxed)

for(i in 1:length(Models_FullRelaxed)){
  k.rate <- max(Models_FullRelaxed[[i]]$index.mat, na.rm = TRUE)
  AICc <- round(Models_FullRelaxed[[i]]$AICc, 2)
  MeanRate <- round(mean(Models_FullRelaxed[[i]]$solution, na.rm = TRUE),2)
  AllTable_FullRelaxed$k.rate[i] <- print(k.rate)
  AllTable_FullRelaxed$AICc[i] <- print(AICc)
  AllTable_FullRelaxed_AICw <- aicw(AllTable_FullRelaxed$AICc)
  AllTable_FullRelaxed$AICw <- round(AllTable_FullRelaxed_AICw$w,2)
  AllTable_FullRelaxed$LogL <- All_Relaxed_Full$LogL_Relaxed_Full
  AllTable_FullRelaxed$MeanRate[i] <- print(MeanRate)
  View(AllTable_FullRelaxed)
}

## Strict Full table
AICc_Strict_Full <- c(Strict_ER_Dual_NoHidden_Model$AICc, 
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

LogL_Strict_Full <- c(Strict_ER_Dual_NoHidden_Model$loglik, 
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


All_Strict_AICc_Full <-as.data.frame(AICc_Strict_Full, row.names = c("StrictFull_ER_Dual_NoHidden_Model", 
                                                                     "StrictFull_ER_NoDual_NoHidden_Model",
                                                                     "StrictFull_ER_DualAllowed_HiddenRatCat2_Model", 
                                                                     "StrictFull_ER_NoDual_HiddenRatCat2_Model", 
                                                                     "StrictFull_ARD_DualAllowed_HiddenRatCat2_Model", 
                                                                     "StrictFull_ARD_Dual_NoHidden_Model", 
                                                                     "StrictFull_ARD_NoDual_HiddenRatCat2_Model", 
                                                                     "StrictFull_ARD_NoDual_NoHidden_Model",
                                                                     "StrictFull_SYM_DualAllowed_HiddenRatCat2_Model",
                                                                     "StrictFull_SYM_Dual_NoHidden_Model",
                                                                     "StrictFull_SYM_NoDual_HiddenRatCat2_Model",
                                                                     "StrictFull_SYM_NoDual_NoHidden_Model"))

All_Strict_Full <-as.data.frame(LogL_Strict_Full, row.names = c("StrictFull_ER_Dual_NoHidden_Model", 
                                                                "StrictFull_ER_NoDual_NoHidden_Model",
                                                                "StrictFull_ER_DualAllowed_HiddenRatCat2_Model", 
                                                                "StrictFull_ER_NoDual_HiddenRatCat2_Model", 
                                                                "StrictFull_ARD_DualAllowed_HiddenRatCat2_Model", 
                                                                "StrictFull_ARD_Dual_NoHidden_Model", 
                                                                "StrictFull_ARD_NoDual_HiddenRatCat2_Model", 
                                                                "StrictFull_ARD_NoDual_NoHidden_Model",
                                                                "StrictFull_SYM_DualAllowed_HiddenRatCat2_Model",
                                                                "StrictFull_SYM_Dual_NoHidden_Model",
                                                                "StrictFull_SYM_NoDual_HiddenRatCat2_Model",
                                                                "StrictFull_SYM_NoDual_NoHidden_Model"))


Models_FullStrict <- c(list(Strict_ARD_Dual_NoHidden_Model,
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

AllTable_FullStrict <- matrix(0, length(Models_FullStrict), 4)
rownames(AllTable_FullStrict) <- rep(c("Full_Strict_ARD_Dual_NoHidden_Model",
                                       "Full_Strict_ARD_DualAllowed_HiddenRatCat2_Model",
                                       "Full_Strict_ARD_NoDual_HiddenRatCat2_Model",
                                       "Full_Strict_ARD_NoDual_NoHidden_Model",
                                       "Full_Strict_ER_Dual_NoHidden_Model",
                                       "Full_Strict_ER_DualAllowed_HiddenRatCat2_Model",
                                       "Full_Strict_ER_NoDual_HiddenRatCat2_Model",
                                       "Full_Strict_ER_NoDual_NoHidden_Model",
                                       "Full_Strict_SYM_Dual_NoHidden_Model",
                                       "Full_Strict_SYM_DualAllowed_HiddenRatCat2_Model",
                                       "Full_Strict_SYM_NoDual_HiddenRatCat2_Model",
                                       "Full_Strict_SYM_NoDual_NoHidden_Model"))

colnames(AllTable_FullStrict) <- c("k.rate", "AICc", "MeanRate", "LogL")

AllTable_FullStrict_Print <- as.data.frame(AllTable_FullStrict)

#Nice table
for(i in 1:length(Models_FullStrict)){
  k.rate <- max(Models_FullStrict[[i]]$index.mat, na.rm = TRUE)
  AICc <- round(Models_FullStrict[[i]]$AICc, 2)
  MeanRate <- round(mean(Models_FullStrict[[i]]$solution, na.rm = TRUE),2)
  AllTable_FullStrict_Print$k.rate[i] <- print(k.rate)
  AllTable_FullStrict_Print$AICc[i] <- print(AICc)
  AllTable_FullStrict_Print_AICw <- aicw(AllTable_FullStrict_Print$AICc)
  AllTable_FullStrict_Print$AICw <- round(AllTable_FullStrict_Print_AICw$w,2)
  AllTable_FullStrict_Print$LogL <- All_Strict_Full$LogL_Strict_Full
  AllTable_FullStrict_Print$MeanRate[i] <- print(MeanRate)
  View(AllTable_FullStrict_Print)
}

