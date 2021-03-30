#SSE models for Strict Full datasets

rm(list=ls())

library(hisse)
library(diversitree)
library(geiger)
library(phytools)
library(ape)
library(tidyverse)

## Read in trait datasets
Full.Strict.Dataset<-read.csv("FullStrict_TraitAnalysis.csv")
Full.Relaxed.Dataset<-read.csv("Full.Relaxed.Dataset.csv")

## Read in phylogeny
FullTree <- read.nexus(file = "FullTree.tre")

## Sampling fractions
strict.arboreal.sampling.fraction <- c(0.48, 0.48, 0.40, 0.50)
relaxed.arboreal.sampling.fraction <- c(0.47, 0.43, 0.47, 0.52)

## Strict with full dataset

## MuSSE null -- no dual transitions
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
Strict.Null.MuS <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1,1,1,1),
                     eps=c(1,1,1,1), root.p = c(1,0,0,0),hidden.states=FALSE,
                     trans.rate=trans.rate)

## MuSSE true-- no dual transitions NEED TO RUN THIS ONE AGAIN
Strict.True.MuS <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1,2,3,4),
                                        eps=c(1,1,1,1), root.p = c(1,0,0,0), hidden.states=FALSE,
                                        trans.rate=trans.rate)

## Character-dependent MuHiSSE model with two hidden states
trans.rate.CharDep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1)
Strict.MuH.CD2 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1,2,3,4,5,6,7,8),
                   eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                   trans.rate=trans.rate.CharDep.TwoHidden)

## Character-independent MuHiSSE model with two hidden states 
trans.rate.CharIndep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
Strict.MuH.CID2 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1,1,1,1,2,2,2,2),
                   eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                   trans.rate=trans.rate.CharIndep.TwoHidden)

## Character-dependent MuHiSSE model with three hidden states
trans.rate.CharDep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2)
Strict.MuH.CD3 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1,2,3,4,5,6,7,8,9,10,11,12),
                                                  eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                  trans.rate=trans.rate.CharDep.ThreeHidden)

## Character-independent MuHiSSE model with three hidden states 
trans.rate.CharIndep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2, make.null=TRUE)
Strict.MuH.CID3 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4)),
                                                    eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharIndep.ThreeHidden)

## Character-dependent MuHiSSE model with four hidden states
trans.rate.CharDep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3)
Strict.MuH.CD4 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1:16),
                                                    eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharDep.FourHidden)

## Character-independent MuHiSSE model with four hidden states 
trans.rate.CharIndep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3, make.null=TRUE)
Strict.MuH.CID4 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4)),
                                                      eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                      trans.rate=trans.rate.CharIndep.FourHidden)

## Character-dependent MuHiSSE model with five hidden states
trans.rate.CharDep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4)
Strict.MuH.CD5 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1:20),
                                                   eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                   trans.rate=trans.rate.CharDep.FiveHidden)

## Character-independent MuHiSSE model with five hidden states 
trans.rate.CharIndep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4, make.null=TRUE)
Strict.MuH.CID5 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4)),
                                                     eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                     trans.rate=trans.rate.CharIndep.FiveHidden)

## Character-dependent MuHiSSE model with six hidden states
trans.rate.CharDep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5)
Strict.MuH.CD6 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1:24),
                                                   eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                   trans.rate=trans.rate.CharDep.SixHidden)

## Character-independent MuHiSSE model with six hidden states 
trans.rate.CharIndep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5, make.null=TRUE)
Strict.MuH.CID6 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4)),
                                                     eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                     trans.rate=trans.rate.CharIndep.SixHidden)

## Character-dependent MuHiSSE model with seven hidden states
trans.rate.CharDep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6)
Strict.MuH.CD7 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1:28),
                                                  eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                  trans.rate=trans.rate.CharDep.SevenHidden)

## Character-independent MuHiSSE model with seven hidden states 
trans.rate.CharIndep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6, make.null=TRUE)
Strict.MuH.CID7 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7, 4)),
                                                    eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0),hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharIndep.SevenHidden)

## Character-dependent MuHiSSE model with eight hidden states
trans.rate.CharDep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7)
Strict.MuH.CD8 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(1:32),
                                                    eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharDep.EightHidden)

## Character-independent MuHiSSE model with eight hidden states 
trans.rate.CharIndep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7, make.null=TRUE)
Strict.MuH.CID8 <- MuHiSSE(phy=FullTree, data=Full.Strict.Dataset, f=strict.arboreal.sampling.fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4),rep(7,4), rep(8,4)),
                                                      eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                      trans.rate=trans.rate.CharIndep.EightHidden)

## Examine model scores

AICc_Strict_SSE_Full <- c(Strict.MuH.CD2$AICc,  Strict.MuH.CD3$AICc ,   
                           Strict.MuH.CD4$AICc,  Strict.MuH.CD5$AICc,    
                           Strict.MuH.CD6$AICc,  Strict.MuH.CD7$AICc,    
                           Strict.MuH.CD8$AICc,  Strict.MuH.CID2$AICc,   
                           Strict.MuH.CID3$AICc, Strict.MuH.CID4$AICc,   
                           Strict.MuH.CID5$AICc ,Strict.MuH.CID6$AICc,   
                           Strict.MuH.CID7$AICc, Strict.MuH.CID8$AICc,   
                           Strict.Null.MuS$AICc, Strict.True.MuS$AICc)

LogL_Strict_SSE_Full <- c(Strict.MuH.CD2$loglik,  Strict.MuH.CD3$loglik ,   
                           Strict.MuH.CD4$loglik,  Strict.MuH.CD5$loglik,    
                           Strict.MuH.CD6$loglik,  Strict.MuH.CD7$loglik,    
                           Strict.MuH.CD8$loglik,  Strict.MuH.CID2$loglik,   
                           Strict.MuH.CID3$loglik, Strict.MuH.CID4$loglik,   
                           Strict.MuH.CID5$loglik ,Strict.MuH.CID6$loglik,   
                           Strict.MuH.CID7$loglik, Strict.MuH.CID8$loglik,   
                           Strict.Null.MuS$loglik, Strict.True.MuS$loglik)


Res_Strict_AICc_SSE_Full <-as.data.frame(AICc_Strict_SSE_Full, row.names = c("Strict.MuH.CD2",  "Strict.MuH.CD3" ,   
                                                                               "Strict.MuH.CD4",  "Strict.MuH.CD5",    
                                                                               "Strict.MuH.CD6",  "Strict.MuH.CD7",    
                                                                               "Strict.MuH.CD8",  "Strict.MuH.CID2",   
                                                                               "Strict.MuH.CID3", "Strict.MuH.CID4",   
                                                                               "Strict.MuH.CID5" ,"Strict.MuH.CID6",   
                                                                               "Strict.MuH.CID7", "Strict.MuH.CID8",   
                                                                               "Strict.Null.MuS", "Strict.True.MuS"))

Res_Strict_LogL_SSE_Full <-as.data.frame(LogL_Strict_SSE_Full, row.names = c("Strict.MuH.CD2",  "Strict.MuH.CD3" ,   
                                                                               "Strict.MuH.CD4",  "Strict.MuH.CD5",    
                                                                               "Strict.MuH.CD6",  "Strict.MuH.CD7",    
                                                                               "Strict.MuH.CD8",  "Strict.MuH.CID2",   
                                                                               "Strict.MuH.CID3", "Strict.MuH.CID4",   
                                                                               "Strict.MuH.CID5" ,"Strict.MuH.CID6",   
                                                                               "Strict.MuH.CID7", "Strict.MuH.CID8",   
                                                                               "Strict.Null.MuS", "Strict.True.MuS"))


Res_Strict_LogL_SSE_Full$AICc<-Res_Strict_AICc_SSE_Full #Strict_SYM_NoDual_HiddenRatCat2_Model wins at 98.2% wAIC
Res_Strict_LogL_SSE_Full$AIC <- Res_Strict_LogL_SSE_Full$AICc$AICc_Strict_SSE_Full
Res_Strict_LogL_SSE_Full$AICw <- aicw(AICc_Strict_SSE_Full)

StrictFullModels <- c(Strict.MuH.CD2,  Strict.MuH.CD3 ,   
                          Strict.MuH.CD4,  Strict.MuH.CD5,    
                          Strict.MuH.CD6,  Strict.MuH.CD7,    
                          Strict.MuH.CD8,  Strict.MuH.CID2,   
                          Strict.MuH.CID3, Strict.MuH.CID4,   
                          Strict.MuH.CID5 ,Strict.MuH.CID6,   
                          Strict.MuH.CID7, Strict.MuH.CID8,   
                          Strict.Null.MuS, Strict.True.MuS)

StrictNames <- c("Strict.MuH.CD2",  "Strict.MuH.CD3" ,   
                 "Strict.MuH.CD4",  "Strict.MuH.CD5",    
                 "Strict.MuH.CD6",  "Strict.MuH.CD7",    
                 "Strict.MuH.CD8",  "Strict.MuH.CID2",   
                 "Strict.MuH.CID3", "Strict.MuH.CID4",   
                 "Strict.MuH.CID5" ,"Strict.MuH.CID6",   
                 "Strict.MuH.CID7", "Strict.MuH.CID8",   
                 "Strict.Null.MuS", "Strict.True.MuS")

#Best fit model for Full Strict dataset is Strict.MuH.CID7
## Adaptive sampling
Strict.Full.MuH.CID7.AdaptiveSampling <- SupportRegionMuHiSSE(muhisse.obj = Strict.MuH.CID7, n.points = 10000)

## Extracting transition rate estimates

StrictFullSolution <- Strict.MuH.CID7$solution

## Transpose the dataframe
StrictFullSolution.Transposed <- as.data.frame(t(as.matrix(StrictFullSolution)))

## a little help from dplyr to keep only the columns with underscores (transition rate estimates)
StrictFullSolution.Rates <- StrictFullSolution.Transposed %>% select(contains("_"))

## 1
StrictFullSolution.Rates1 <- StrictFullSolution.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
StrictFullSolution.Rates.00.10 <- StrictFullSolution.Rates1 %>% select(contains("q00"))
## q10 > q00
StrictFullSolution.Rates.10.00 <- StrictFullSolution.Rates1 %>% select(contains("q10"))

## 2
StrictFullSolution.Rates.q01q00 <- StrictFullSolution.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
StrictFullSolution.Rates.00.01 <- StrictFullSolution.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
StrictFullSolution.Rates.01.00 <- StrictFullSolution.Rates.q01q00 %>% select(contains("q01"))

## 3
StrictFullSolution.Rates.q01q11 <- StrictFullSolution.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
StrictFullSolution.Rates.01.11 <- StrictFullSolution.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
StrictFullSolution.Rates.11.01 <- StrictFullSolution.Rates.q01q11 %>% select(contains("q11"))

## 4
StrictFullSolution.Rates.q10q11 <- StrictFullSolution.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
StrictFullSolution.Rates.10.11 <- StrictFullSolution.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
StrictFullSolution.Rates.11.10 <- StrictFullSolution.Rates.q10q11 %>% select(contains("q11"))

## Combine all into a new dataframe for plotting
StrictFull_qrates <- list(StrictFullSolution.Rates.00.10, 
            StrictFullSolution.Rates.10.00, 
            StrictFullSolution.Rates.00.01, 
            StrictFullSolution.Rates.01.00, 
            StrictFullSolution.Rates.01.11, 
            StrictFullSolution.Rates.11.01, 
            StrictFullSolution.Rates.10.11, 
            StrictFullSolution.Rates.11.10) %>%
  setNames(c("00.10",
             "10.00",
             "00.01",
             "01.00",
             "01.11",
             "11.01",
             "10.11",
             "11.10"))



#Just sample one
for(i in 1:length(StrictFull_qrates)){
  MuHiSSEStrictFullRates <- sample(StrictFull_qrates[[i]][1]) * 100
  MuHiSSEStrictFullRatesPrint <- print(round(MuHiSSEStrictFullRates, 4))
  #write.csv(MuHiSSEStrictFullRates[[i]][1], "MuHiSSEStrictFullRates.Strict.MuH.CID7.csv")
}

## Ancestral State Reconstruction
## Reconstruct ancestral states for this model

Strict.Full.MuH.CID7.recon <-
  hisse::MarginReconMuHiSSE(
    phy = Strict.MuH.CID7$phy,
    data = Strict.MuH.CID7$data,
    f = Strict.MuH.CID7$f,
    pars = Strict.MuH.CID7$solution,
    hidden.states = 7,
    aic = Strict.MuH.CID7$AICc,
    n.cores = 3,
    verbose = TRUE
  )

## PLOT
Strict.Full.MuH.CID7.recon 
Strict.Full.MuH.CID7.proc <- m_process_recon(Strict.Full.MuH.CID7.recon)
Strict.Full.MuH.CID7.recon.AncestralStates <- m_trait_recon(
processed_recon = Strict.Full.MuH.CID7.proc,
states_of_first_character = c('0','1'),
states_of_second_character = c('0','1'),
cutoff = as.numeric(c('0.5','0.5')),
colors = c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"),
show_tip_labels = FALSE,
tree_layout = 'fan',
tree_direction = 'up',
time_axis_ticks = 10,
open_angle = 10)

## Plot with clade labels for clades with padbearing + aboreal taxa
Strict.Full.MuH.CID7.recon.AncestralStates + geom_cladelabel(node=2701, label="Gekkota", align=T, geom='label', fill='white',offset.text = 6)+
  geom_cladelabel(node=3426, label="Scincidae", align=T, geom='label', fill='white',offset.text = 40)+
  geom_cladelabel(node=5190, label="Dactyloidae", align=T, geom='label', fill='white',offset.text = 6)

## Net Diverisification
StrictFull.MuH.CID7.NetDiv <- m_scatterplot(
  processed_recon = Strict.Full.MuH.CID7.proc,
  parameter = 'net.div',
  states_names = c('00','01','10','11'), colors = c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"), 
  plot_as_waiting_time = FALSE) + theme_classic()

StrictFull.MuH.CID7.NetDiv + ylab("Net Diversfication")

