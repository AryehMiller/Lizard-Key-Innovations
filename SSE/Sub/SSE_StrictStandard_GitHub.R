## SSE models for the Strict Standard dataset

rm(list=ls())

library(hisse)
library(diversitree)
library(geiger)
library(phytools)
library(ape)
library(tidyverse)
library(geiger)

## Read in trait datasets (strict standard)
Standard.Strict.Dataset<-read.csv("StandardStrict_TraitAnalysis.csv")

## Read in phylogeny
StandardTree <- read.nexus(file = "StandardTree.tre")

### Sampling fractions
Standard.Strict.Sampling.Fraction <- c(0.40, 0.48, 0.40, 0.50)

## MuSSE null -- no dual transitions
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
Strict.Standard.Null.MuS <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1,1,1,1),
                     eps=c(1,1,1,1), root.p = c(1,0,0,0),hidden.states=FALSE,
                     trans.rate=trans.rate)

## MuSSE true-- no dual transitions NEED TO RUN THIS ONE AGAIN
Strict.Standard.True.MuS <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1,2,3,4),
                                        eps=c(1,1,1,1), root.p = c(1,0,0,0), hidden.states=FALSE,
                                        trans.rate=trans.rate)


## Character-dependent MuHiSSE model with two hidden states
trans.rate.CharDep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1)
Strict.Standard.MuH.CD2 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1,2,3,4,5,6,7,8),
                   eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                   trans.rate=trans.rate.CharDep.TwoHidden)

## Character-independent MuHiSSE model with two hidden states 
trans.rate.CharIndep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
Strict.Standard.MuH.CID2 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1,1,1,1,2,2,2,2),
                   eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                   trans.rate=trans.rate.CharIndep.TwoHidden)

## Character-dependent MuHiSSE model with three hidden states
trans.rate.CharDep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2)
Strict.Standard.MuH.CD3 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1,2,3,4,5,6,7,8,9,10,11,12),
                                                  eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                  trans.rate=trans.rate.CharDep.ThreeHidden)

## Character-independent MuHiSSE model with three hidden states 
trans.rate.CharIndep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2, make.null=TRUE)
Strict.Standard.MuH.CID3 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4)),
                                                    eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharIndep.ThreeHidden)

## Character-dependent MuHiSSE model with four hidden states
trans.rate.CharDep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3)
Strict.Standard.MuH.CD4 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1:16),
                                                    eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharDep.FourHidden)

## Character-independent MuHiSSE model with four hidden states 
trans.rate.CharIndep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3, make.null=TRUE)
Strict.Standard.MuH.CID4 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4)),
                                                      eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                      trans.rate=trans.rate.CharIndep.FourHidden)

## Character-dependent MuHiSSE model with five hidden states
trans.rate.CharDep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4)
Strict.Standard.MuH.CD5 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1:20),
                                                   eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                   trans.rate=trans.rate.CharDep.FiveHidden)

## Character-independent MuHiSSE model with five hidden states 
trans.rate.CharIndep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4, make.null=TRUE)
Strict.Standard.MuH.CID5 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4)),
                                                     eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                     trans.rate=trans.rate.CharIndep.FiveHidden)

## Character-dependent MuHiSSE model with six hidden states
trans.rate.CharDep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5)
Strict.Standard.MuH.CD6 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1:24),
                                                   eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                   trans.rate=trans.rate.CharDep.SixHidden)

## Character-independent MuHiSSE model with six hidden states 
trans.rate.CharIndep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5, make.null=TRUE)
Strict.Standard.MuH.CID6 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4)),
                                                     eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                     trans.rate=trans.rate.CharIndep.SixHidden)

## Character-dependent MuHiSSE model with seven hidden states
trans.rate.CharDep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6)
Strict.Standard.MuH.CD7 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1:28),
                                                  eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                  trans.rate=trans.rate.CharDep.SevenHidden)

## Character-independent MuHiSSE model with seven hidden states 
trans.rate.CharIndep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6, make.null=TRUE)
Strict.Standard.MuH.CID7 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7, 4)),
                                                    eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0),hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharIndep.SevenHidden)

## Character-dependent MuHiSSE model with eight hidden states
trans.rate.CharDep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7)
Strict.Standard.MuH.CD8 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(1:32),
                                                    eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharDep.EightHidden)

## Character-independent MuHiSSE model with eight hidden states 
trans.rate.CharIndep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7, make.null=TRUE)
Strict.Standard.MuH.CID8 <- MuHiSSE(phy=StandardTree, data=Standard.Strict.Dataset, f=Standard.Strict.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4),rep(7,4), rep(8,4)),
                                                      eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                      trans.rate=trans.rate.CharIndep.EightHidden)

## Examine model scores

AICc_Strict_Standard_SSE_Full <- c(Strict.Standard.MuH.CD2$AICc,  Strict.Standard.MuH.CD3$AICc ,   
                                   Strict.Standard.MuH.CD4$AICc,  Strict.Standard.MuH.CD5$AICc,    
                                   Strict.Standard.MuH.CD6$AICc,  Strict.Standard.MuH.CD7$AICc,    
                                   Strict.Standard.MuH.CD8$AICc,  Strict.Standard.MuH.CID2$AICc,   
                                   Strict.Standard.MuH.CID3$AICc, Strict.Standard.MuH.CID4$AICc,   
                                   Strict.Standard.MuH.CID5$AICc ,Strict.Standard.MuH.CID6$AICc,   
                                   Strict.Standard.MuH.CID7$AICc, Strict.Standard.MuH.CID8$AICc,   
                                   Strict.Standard.Null.MuS$AICc, Strict.Standard.True.MuS$AICc)

LogL_Strict_Standard_SSE_Full <- c(Strict.Standard.MuH.CD2$loglik,  Strict.Standard.MuH.CD3$loglik ,   
                                   Strict.Standard.MuH.CD4$loglik,  Strict.Standard.MuH.CD5$loglik,    
                                   Strict.Standard.MuH.CD6$loglik,  Strict.Standard.MuH.CD7$loglik,    
                                   Strict.Standard.MuH.CD8$loglik,  Strict.Standard.MuH.CID2$loglik,   
                                   Strict.Standard.MuH.CID3$loglik, Strict.Standard.MuH.CID4$loglik,   
                                   Strict.Standard.MuH.CID5$loglik ,Strict.Standard.MuH.CID6$loglik,   
                                   Strict.Standard.MuH.CID7$loglik, Strict.Standard.MuH.CID8$loglik,   
                                   Strict.Standard.Null.MuS$loglik, Strict.Standard.True.MuS$loglik)


Res_Strict_Standard_AICc_SSE_Full <-as.data.frame(AICc_Strict_Standard_SSE_Full, row.names = c("Strict.Standard.MuH.CD2",  "Strict.Standard.MuH.CD3" ,   
                                                                                               "Strict.Standard.MuH.CD4",  "Strict.Standard.MuH.CD5",    
                                                                                               "Strict.Standard.MuH.CD6",  "Strict.Standard.MuH.CD7",    
                                                                                               "Strict.Standard.MuH.CD8",  "Strict.Standard.MuH.CID2",   
                                                                                               "Strict.Standard.MuH.CID3", "Strict.Standard.MuH.CID4",   
                                                                                               "Strict.Standard.MuH.CID5" ,"Strict.Standard.MuH.CID6",   
                                                                                               "Strict.Standard.MuH.CID7", "Strict.Standard.MuH.CID8",   
                                                                                               "Strict.Standard.Null.MuS", "Strict.Standard.True.MuS"))

Res_Strict_Standard_LogL_SSE_Full <-as.data.frame(LogL_Strict_Standard_SSE_Full, row.names = c("Strict.Standard.MuH.CD2",  "Strict.Standard.MuH.CD3" ,   
                                                                                               "Strict.Standard.MuH.CD4",  "Strict.Standard.MuH.CD5",    
                                                                                               "Strict.Standard.MuH.CD6",  "Strict.Standard.MuH.CD7",    
                                                                                               "Strict.Standard.MuH.CD8",  "Strict.Standard.MuH.CID2",   
                                                                                               "Strict.Standard.MuH.CID3", "Strict.Standard.MuH.CID4",   
                                                                                               "Strict.Standard.MuH.CID5" ,"Strict.Standard.MuH.CID6",   
                                                                                               "Strict.Standard.MuH.CID7", "Strict.Standard.MuH.CID8",   
                                                                                               "Strict.Standard.Null.MuS", "Strict.Standard.True.MuS"))


Res_Strict_Standard_LogL_SSE_Full$AICc<-Res_Strict_Standard_AICc_SSE_Full #Strict_SYM_NoDual_HiddenRatCat2_Model wins at 98.2% wAIC
Res_Strict_Standard_LogL_SSE_Full$AIC <- Res_Strict_Standard_LogL_SSE_Full$AICc$AICc_Strict_Standard_SSE_Full
Res_Strict_Standard_LogL_SSE_Full$AICw <- aicw(AICc_Strict_Standard_SSE_Full)

##write.csv here if you want

## Best fit model for strict dataset is Strict.Standard.MuH.CID5
## Reconstruct ancestral states for this model

Strict.Standard.MuH.CID5.recon <-
  hisse::MarginReconMuHiSSE(
    phy = Strict.Standard.MuH.CID5$phy,
    data = Strict.Standard.MuH.CID5$data,
    f = Strict.Standard.MuH.CID5$f,
    pars = Strict.Standard.MuH.CID5$solution,
    hidden.states = 5,
    aic = Strict.Standard.MuH.CID5$AICc,
    n.cores = 3,
    verbose = TRUE
  )

## PLOTTING
## Ancestral State Reconstruction
Strict.Standard.MuH.CID5.recon 
Strict.Standard.MuH.CID5.proc <- m_process_recon(Strict.Standard.MuH.CID5.recon)
Strict.Standard.MuH.CID5.recon.AncestralStates <- m_trait_recon(
  processed_recon = Strict.Standard.MuH.CID5.proc,
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
Strict.Standard.MuH.CID5.recon.AncestralStates + geom_cladelabel(node=2362, label="Gekkota", align=T, geom='label', fill='white',offset.text = 6)+
  geom_cladelabel(node=3052, label="Scincidae", align=T, geom='label', fill='white',offset.text = 40)+
  geom_cladelabel(node=4526, label="Dactyloidae", align=T, geom='label', fill='white',offset.text = 6)

## Net Diverisification
StrictFull.MuH.CID5.NetDiv <- m_scatterplot(
  processed_recon = Strict.Standard.MuH.CID5.proc,
  parameter = 'net.div',
  states_names = c('00','01','10','11'), colors = c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"), 
  plot_as_waiting_time = FALSE) + theme_classic()

StrictFull.MuH.CID5.NetDiv + ylab("Net Diversfication")

## Transition rates
Strict.Standard.MuH.CID5

StrictStandardSolution <- Strict.Standard.MuH.CID5$solution

## Transpose the dataframe
StrictStandardSolution.Transposed <- as.data.frame(t(as.matrix(StrictStandardSolution)))

## a little help from dplyr to keep only the columns with underscores (transition rate estimates)
StrictStandardSolution.Rates <- StrictStandardSolution.Transposed %>% select(contains("_"))

## 1
StrictStandardSolution.Rates1 <- StrictStandardSolution.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
StrictStandardSolution.Rates.00.10 <- StrictStandardSolution.Rates1 %>% select(contains("q00"))
## q10 > q00
StrictStandardSolution.Rates.10.00 <- StrictStandardSolution.Rates1 %>% select(contains("q10"))

## 2
StrictStandardSolution.Rates.q01q00 <- StrictStandardSolution.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
StrictStandardSolution.Rates.00.01 <- StrictStandardSolution.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
StrictStandardSolution.Rates.01.00 <- StrictStandardSolution.Rates.q01q00 %>% select(contains("q01"))

## 3
StrictStandardSolution.Rates.q01q11 <- StrictStandardSolution.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
StrictStandardSolution.Rates.01.11 <- StrictStandardSolution.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
StrictStandardSolution.Rates.11.01 <- StrictStandardSolution.Rates.q01q11 %>% select(contains("q11"))

## 4
StrictStandardSolution.Rates.q10q11 <- StrictStandardSolution.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
StrictStandardSolution.Rates.10.11 <- StrictStandardSolution.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
StrictStandardSolution.Rates.11.10 <- StrictStandardSolution.Rates.q10q11 %>% select(contains("q11"))

## Collate all into a new dataframe for plotting
StrictStandard_qrates <- list(StrictStandardSolution.Rates.00.10, 
                          StrictStandardSolution.Rates.10.00, 
                          StrictStandardSolution.Rates.00.01, 
                          StrictStandardSolution.Rates.01.00, 
                          StrictStandardSolution.Rates.01.11, 
                          StrictStandardSolution.Rates.11.01, 
                          StrictStandardSolution.Rates.10.11, 
                          StrictStandardSolution.Rates.11.10) %>%
  setNames(c("00.10",
             "10.00",
             "00.01",
             "01.00",
             "01.11",
             "11.01",
             "10.11",
             "11.10"))



#Just sample one rate for each transition category (so all rates are the same across hidden states so this is OK)
for(i in 1:length(StrictStandard_qrates)){
  MuHiSSEStrictStandardRates <- sample(StrictStandard_qrates[[i]][1]) * 100
  MuHiSSEStrictStandardRatesPrint <- print(round(MuHiSSEStrictStandardRates, 4))
  #write.csv(MuHiSSEStrictStandardRates[[i]][1], "MuHiSSEStrictStandardRates.Strict.MuH.CID7.csv")
}
