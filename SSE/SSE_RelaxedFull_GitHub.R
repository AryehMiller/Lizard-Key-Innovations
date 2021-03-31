#SSE models Relaxed Full

rm(list=ls())

library(hisse)
library(diversitree)
library(geiger)
library(phytools)
library(ape)
library(tidyverse)

## Read in trait dataset

Full.Relaxed.Dataset<-read.csv("FullRelaxed_TraitAnalysis.csv")

## Read in phylogeny
FullTree <- read.nexus(file = "FullTree.tre")

## Sampling fractions
Relaxed.Full.Sampling.Fraction <- c(0.47, 0.43, 0.47, 0.52)

## MuSSE null -- no dual transitions
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
Relaxed.Null.MuS <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1,1,1,1),
                           eps=c(1,1,1,1), root.p = c(1,0,0,0),hidden.states=FALSE,
                           trans.rate=trans.rate)

## MuSSE true-- no dual transitions
Relaxed.True.MuS <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1,2,3,4),
                           eps=c(1,1,1,1), root.p = c(1,0,0,0), hidden.states=FALSE,
                           trans.rate=trans.rate)

## Character-dependent MuHiSSE model with two hidden states
trans.rate.CharDep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1)
Relaxed.MuH.CD2 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1,2,3,4,5,6,7,8),
                          eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                          trans.rate=trans.rate.CharDep.TwoHidden)

## Character-independent MuHiSSE model with two hidden states 
trans.rate.CharIndep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
Relaxed.MuH.CID2 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1,1,1,1,2,2,2,2),
                           eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                           trans.rate=trans.rate.CharIndep.TwoHidden)

## Character-dependent MuHiSSE model with three hidden states
trans.rate.CharDep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2)
Relaxed.MuH.CD3 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1,2,3,4,5,6,7,8,9,10,11,12),
                          eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                          trans.rate=trans.rate.CharDep.ThreeHidden)

## Character-independent MuHiSSE model with three hidden states 
trans.rate.CharIndep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2, make.null=TRUE)
Relaxed.MuH.CID3 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4)),
                           eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                           trans.rate=trans.rate.CharIndep.ThreeHidden)

## Character-dependent MuHiSSE model with four hidden states
trans.rate.CharDep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3)
Relaxed.MuH.CD4 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1:16),
                          eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                          trans.rate=trans.rate.CharDep.FourHidden)

## Character-independent MuHiSSE model with four hidden states 
trans.rate.CharIndep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3, make.null=TRUE)
Relaxed.MuH.CID4 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4)),
                           eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                           trans.rate=trans.rate.CharIndep.FourHidden)

save.image(file = "SSE_Relaxed.RData")

## Character-dependent MuHiSSE model with five hidden states
trans.rate.CharDep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4)
Relaxed.MuH.CD5 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1:20),
                          eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                          trans.rate=trans.rate.CharDep.FiveHidden)

## Character-independent MuHiSSE model with five hidden states 
trans.rate.CharIndep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4, make.null=TRUE)
Relaxed.MuH.CID5 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4)),
                           eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                           trans.rate=trans.rate.CharIndep.FiveHidden)

## Character-dependent MuHiSSE model with six hidden states
trans.rate.CharDep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5)
Relaxed.MuH.CD6 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1:24),
                          eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                          trans.rate=trans.rate.CharDep.SixHidden)

## Character-independent MuHiSSE model with six hidden states 
trans.rate.CharIndep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5, make.null=TRUE)
Relaxed.MuH.CID6 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4)),
                           eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                           trans.rate=trans.rate.CharIndep.SixHidden)

## Character-dependent MuHiSSE model with seven hidden states
trans.rate.CharDep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6)
Relaxed.MuH.CD7 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1:28),
                          eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                          trans.rate=trans.rate.CharDep.SevenHidden)

## Character-independent MuHiSSE model with seven hidden states 
trans.rate.CharIndep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6, make.null=TRUE)
Relaxed.MuH.CID7 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7, 4)),
                           eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0),hidden.states=TRUE,
                           trans.rate=trans.rate.CharIndep.SevenHidden)

save.image(file = "SSE_Relaxed.RData")


## Character-dependent MuHiSSE model with eight hidden states
trans.rate.CharDep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7)
Relaxed.MuH.CD8 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(1:32),
                          eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                          trans.rate=trans.rate.CharDep.EightHidden)

## Character-independent MuHiSSE model with eight hidden states 
trans.rate.CharIndep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7, make.null=TRUE)
Relaxed.MuH.CID8 <- MuHiSSE(phy=FullTree, data=Full.Relaxed.Dataset, f=Relaxed.Full.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4),rep(7,4), rep(8,4)),
                           eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                           trans.rate=trans.rate.CharIndep.EightHidden)

## Relaxed Full Analysis

## State-Averaged Transition Rates 
RelaxedFullSolution <- Relaxed.MuH.CID8$solution

## Transpose the dataframe
RelaxedFullSolution.Transposed <- as.data.frame(t(as.matrix(RelaxedFullSolution)))

## a little help from dplyr to keep only the columns with underscores (transition rate estimates)
RelaxedFullSolution.Rates <- RelaxedFullSolution.Transposed %>% select(contains("_"))

## 1
RelaxedFullSolution.Rates1 <- RelaxedFullSolution.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
RelaxedFullSolution.Rates.00.10 <- RelaxedFullSolution.Rates1 %>% select(contains("q00"))
## q10 > q00
RelaxedFullSolution.Rates.10.00 <- RelaxedFullSolution.Rates1 %>% select(contains("q10"))

## 2
RelaxedFullSolution.Rates.q01q00 <- RelaxedFullSolution.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
RelaxedFullSolution.Rates.00.01 <- RelaxedFullSolution.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
RelaxedFullSolution.Rates.01.00 <- RelaxedFullSolution.Rates.q01q00 %>% select(contains("q01"))

## 3
RelaxedFullSolution.Rates.q01q11 <- RelaxedFullSolution.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
RelaxedFullSolution.Rates.01.11 <- RelaxedFullSolution.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
RelaxedFullSolution.Rates.11.01 <- RelaxedFullSolution.Rates.q01q11 %>% select(contains("q11"))

## 4
RelaxedFullSolution.Rates.q10q11 <- RelaxedFullSolution.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
RelaxedFullSolution.Rates.10.11 <- RelaxedFullSolution.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
RelaxedFullSolution.Rates.11.10 <- RelaxedFullSolution.Rates.q10q11 %>% select(contains("q11"))

## Combine all into a new dataframe for plotting
RelaxedFull.qrates <- list(RelaxedFullSolution.Rates.00.10, 
            RelaxedFullSolution.Rates.10.00, 
            RelaxedFullSolution.Rates.00.01, 
            RelaxedFullSolution.Rates.01.00, 
            RelaxedFullSolution.Rates.01.11, 
            RelaxedFullSolution.Rates.11.01, 
            RelaxedFullSolution.Rates.10.11, 
            RelaxedFullSolution.Rates.11.10) %>%
  setNames(c("00.10",
             "10.00",
             "00.01",
             "01.00",
             "01.11",
             "11.01",
             "10.11",
             "11.10"))



#Just sample one
for(i in 1:length(RelaxedFull.qrates)){
  MuHiSSERelaxedFullRates <- sample(RelaxedFull.qrates[[i]],1)
  print(round(MuHiSSERelaxedFullRates, 4))
}

## Ancestral State Reconstruction

Relaxed.MuH.CID8.recon <-
  hisse::MarginReconMuHiSSE(
    phy = Relaxed.MuH.CID8$phy,
    data = Relaxed.MuH.CID8$data,
    f = Relaxed.MuH.CID8$f,
    pars = Relaxed.MuH.CID8$solution,
    hidden.states = 8,
    aic = Relaxed.MuH.CID8$AICc,
    n.cores = 3,
    verbose = TRUE
  )

m_proc <- m_process_recon(Relaxed.MuH.CID8.recon)
Relaxed.MuH.CID8.recon.AncestralStates <- m_trait_recon(
processed_recon = m_proc,
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
Relaxed.MuH.CID8.recon.AncestralStates + geom_cladelabel(node=2701, label="Gekkota", align=T, geom='label', fill='white',offset.text = 6)+
    geom_cladelabel(node=3426, label="Scincidae", align=T, geom='label', fill='white',offset.text = 40)+
    geom_cladelabel(node=5190, label="Dactyloidae", align=T, geom='label', fill='white',offset.text = 6)
  
  #findMRCA
  ## dev.size("in") 17.95833 10.75000 x 1.3
  ## Net Diverisification
RelaxedScatterFullNetDiv <- m_scatterplot(
    processed_recon = m_proc,
    parameter = 'net.div',
    states_names = c('00','01','10','11'), colors = c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"), 
    plot_as_waiting_time = FALSE) + theme_classic()
  
RelaxedScatterFullNetDiv + ylab("Net Diversfication")

  
  

