#SSE models Standard Relaxed
rm(list=ls())

library(hisse)
library(diversitree)
library(geiger)
library(phytools)
library(ape)
library(tidyverse)

## Read in trait datasets (strict standard)
Standard.Relaxed.Dataset<-read.csv("StandardRelaxed_TraitAnalysis.csv")

## Read in phylogeny
StandardTree <- read.nexus(file = "StandardTree.tre")

## Sampling fraction

Standard.Relaxed.Sampling.Fraction <- c(0.37, 0.43, 0.47, 0.52)

## MuSSE null -- no dual transitions
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
Relaxed.Standard.Null.MuS <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1,1,1,1),
                     eps=c(1,1,1,1), root.p = c(1,0,0,0),hidden.states=FALSE,
                     trans.rate=trans.rate)

## MuSSE true-- no dual transitions NEED TO RUN THIS ONE AGAIN
Relaxed.Standard.True.MuS <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1,2,3,4),
                                        eps=c(1,1,1,1), root.p = c(1,0,0,0), hidden.states=FALSE,
                                        trans.rate=trans.rate)

## Character-dependent MuHiSSE model with two hidden states
trans.rate.CharDep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1)
Relaxed.Standard.MuH.CD2 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1,2,3,4,5,6,7,8),
                   eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                   trans.rate=trans.rate.CharDep.TwoHidden)

## Character-independent MuHiSSE model with two hidden states 
trans.rate.CharIndep.TwoHidden <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
Relaxed.Standard.MuH.CID2 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1,1,1,1,2,2,2,2),
                   eps=rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states=TRUE,
                   trans.rate=trans.rate.CharIndep.TwoHidden)

## Character-dependent MuHiSSE model with three hidden states
trans.rate.CharDep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2)
Relaxed.Standard.MuH.CD3 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1,2,3,4,5,6,7,8,9,10,11,12),
                                                  eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                  trans.rate=trans.rate.CharDep.ThreeHidden)

## Character-independent MuHiSSE model with three hidden states 
trans.rate.CharIndep.ThreeHidden <- TransMatMakerMuHiSSE(hidden.traits=2, make.null=TRUE)
Relaxed.Standard.MuH.CID3 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4)),
                                                    eps=rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharIndep.ThreeHidden)

## Character-dependent MuHiSSE model with four hidden states
trans.rate.CharDep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3)
Relaxed.Standard.MuH.CD4 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1:16),
                                                    eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharDep.FourHidden)

## Character-independent MuHiSSE model with four hidden states 
trans.rate.CharIndep.FourHidden <- TransMatMakerMuHiSSE(hidden.traits=3, make.null=TRUE)
Relaxed.Standard.MuH.CID4 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4)),
                                                      eps=rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                      trans.rate=trans.rate.CharIndep.FourHidden)

## Character-dependent MuHiSSE model with five hidden states
trans.rate.CharDep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4)
Relaxed.Standard.MuH.CD5 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1:20),
                                                   eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                   trans.rate=trans.rate.CharDep.FiveHidden)

## Character-independent MuHiSSE model with five hidden states 
trans.rate.CharIndep.FiveHidden <- TransMatMakerMuHiSSE(hidden.traits=4, make.null=TRUE)
Relaxed.Standard.MuH.CID5 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4)),
                                                     eps=rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                     trans.rate=trans.rate.CharIndep.FiveHidden)

## Character-dependent MuHiSSE model with six hidden states
trans.rate.CharDep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5)
Relaxed.Standard.MuH.CD6 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1:24),
                                                   eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                   trans.rate=trans.rate.CharDep.SixHidden)

## Character-independent MuHiSSE model with six hidden states 
trans.rate.CharIndep.SixHidden <- TransMatMakerMuHiSSE(hidden.traits=5, make.null=TRUE)
Relaxed.Standard.MuH.CID6 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4)),
                                                     eps=rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                     trans.rate=trans.rate.CharIndep.SixHidden)

## Character-dependent MuHiSSE model with seven hidden states
trans.rate.CharDep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6)
Relaxed.Standard.MuH.CD7 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1:28),
                                                  eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                  trans.rate=trans.rate.CharDep.SevenHidden)

## Character-independent MuHiSSE model with seven hidden states 
trans.rate.CharIndep.SevenHidden <- TransMatMakerMuHiSSE(hidden.traits=6, make.null=TRUE)
Relaxed.Standard.MuH.CID7 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7, 4)),
                                                    eps=rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0),hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharIndep.SevenHidden)

## Character-dependent MuHiSSE model with eight hidden states
trans.rate.CharDep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7)
Relaxed.Standard.MuH.CD8 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(1:32),
                                                    eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                    trans.rate=trans.rate.CharDep.EightHidden)

## Character-independent MuHiSSE model with eight hidden states 
trans.rate.CharIndep.EightHidden <- TransMatMakerMuHiSSE(hidden.traits=7, make.null=TRUE)
Relaxed.Standard.MuH.CID8 <- MuHiSSE(phy=StandardTree, data=Standard.Relaxed.Dataset, f=Standard.Relaxed.Sampling.Fraction, turnover=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4),rep(7,4), rep(8,4)),
                                                      eps=rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states=TRUE,
                                                      trans.rate=trans.rate.CharIndep.EightHidden)

## Examine model scores

AICc_Relaxed.Standard_SSE_Full <- c(Relaxed.Standard.MuH.CD2$AICc,  Relaxed.Standard.MuH.CD3$AICc ,   
                           Relaxed.Standard.MuH.CD4$AICc,  Relaxed.Standard.MuH.CD5$AICc,    
                           Relaxed.Standard.MuH.CD6$AICc,  Relaxed.Standard.MuH.CD7$AICc,    
                           Relaxed.Standard.MuH.CD8$AICc,  Relaxed.Standard.MuH.CID2$AICc,   
                           Relaxed.Standard.MuH.CID3$AICc, Relaxed.Standard.MuH.CID4$AICc,   
                           Relaxed.Standard.MuH.CID5$AICc ,Relaxed.Standard.MuH.CID6$AICc,   
                           Relaxed.Standard.MuH.CID7$AICc, Relaxed.Standard.MuH.CID8$AICc,   
                           Relaxed.Standard.Null.MuS$AICc, Relaxed.Standard.True.MuS$AICc)

LogL_Relaxed.Standard_SSE_Full <- c(Relaxed.Standard.MuH.CD2$loglik,  Relaxed.Standard.MuH.CD3$loglik ,   
                           Relaxed.Standard.MuH.CD4$loglik,  Relaxed.Standard.MuH.CD5$loglik,    
                           Relaxed.Standard.MuH.CD6$loglik,  Relaxed.Standard.MuH.CD7$loglik,    
                           Relaxed.Standard.MuH.CD8$loglik,  Relaxed.Standard.MuH.CID2$loglik,   
                           Relaxed.Standard.MuH.CID3$loglik, Relaxed.Standard.MuH.CID4$loglik,   
                           Relaxed.Standard.MuH.CID5$loglik ,Relaxed.Standard.MuH.CID6$loglik,   
                           Relaxed.Standard.MuH.CID7$loglik, Relaxed.Standard.MuH.CID8$loglik,   
                           Relaxed.Standard.Null.MuS$loglik, Relaxed.Standard.True.MuS$loglik)


Res_Relaxed.Standard_AICc_SSE_Full <-as.data.frame(AICc_Relaxed.Standard_SSE_Full, row.names = c("Relaxed.Standard.MuH.CD2",  "Relaxed.Standard.MuH.CD3" ,   
                                                                               "Relaxed.Standard.MuH.CD4",  "Relaxed.Standard.MuH.CD5",    
                                                                               "Relaxed.Standard.MuH.CD6",  "Relaxed.Standard.MuH.CD7",    
                                                                               "Relaxed.Standard.MuH.CD8",  "Relaxed.Standard.MuH.CID2",   
                                                                               "Relaxed.Standard.MuH.CID3", "Relaxed.Standard.MuH.CID4",   
                                                                               "Relaxed.Standard.MuH.CID5" ,"Relaxed.Standard.MuH.CID6",   
                                                                               "Relaxed.Standard.MuH.CID7", "Relaxed.Standard.MuH.CID8",   
                                                                               "Relaxed.Standard.Null.MuS", "Relaxed.Standard.True.MuS"))

Res_Relaxed.Standard_LogL_SSE_Full <-as.data.frame(LogL_Relaxed.Standard_SSE_Full, row.names = c("Relaxed.Standard.MuH.CD2",  "Relaxed.Standard.MuH.CD3" ,   
                                                                               "Relaxed.Standard.MuH.CD4",  "Relaxed.Standard.MuH.CD5",    
                                                                               "Relaxed.Standard.MuH.CD6",  "Relaxed.Standard.MuH.CD7",    
                                                                               "Relaxed.Standard.MuH.CD8",  "Relaxed.Standard.MuH.CID2",   
                                                                               "Relaxed.Standard.MuH.CID3", "Relaxed.Standard.MuH.CID4",   
                                                                               "Relaxed.Standard.MuH.CID5" ,"Relaxed.Standard.MuH.CID6",   
                                                                               "Relaxed.Standard.MuH.CID7", "Relaxed.Standard.MuH.CID8",   
                                                                               "Relaxed.Standard.Null.MuS", "Relaxed.Standard.True.MuS"))


Res_Relaxed.Standard_LogL_SSE_Full$AICc<-Res_Relaxed.Standard_AICc_SSE_Full
Res_Relaxed.Standard_LogL_SSE_Full$AIC <- Res_Relaxed.Standard_LogL_SSE_Full$AICc$AICc_Relaxed_SSE_Standard
Res_Relaxed.Standard_LogL_SSE_Full$AICw <- aicw(AICc_Strict_SSE_Full)
Strict.MuH.CID7 


## Examine model scores

AIC_Relaxed.Standard_SSE_Full <- c(Relaxed.Standard.MuH.CD2$AIC,  Relaxed.Standard.MuH.CD3$AIC ,   
                                   Relaxed.Standard.MuH.CD4$AIC,  Relaxed.Standard.MuH.CD5$AIC,    
                                   Relaxed.Standard.MuH.CD6$AIC,  Relaxed.Standard.MuH.CD7$AIC,    
                                   Relaxed.Standard.MuH.CD8$AIC,  Relaxed.Standard.MuH.CID2$AIC,   
                                   Relaxed.Standard.MuH.CID3$AIC, Relaxed.Standard.MuH.CID4$AIC,   
                                   Relaxed.Standard.MuH.CID5$AIC ,Relaxed.Standard.MuH.CID6$AIC,   
                                   Relaxed.Standard.MuH.CID7$AIC, Relaxed.Standard.MuH.CID8$AIC,   
                                   Relaxed.Standard.Null.MuS$AIC, Relaxed.Standard.True.MuS$AIC)

AICc_Relaxed.Standard_SSE_Full <- c(Relaxed.Standard.MuH.CD2$AICc,  Relaxed.Standard.MuH.CD3$AICc ,   
                                    Relaxed.Standard.MuH.CD4$AICc,  Relaxed.Standard.MuH.CD5$AICc,    
                                    Relaxed.Standard.MuH.CD6$AICc,  Relaxed.Standard.MuH.CD7$AICc,    
                                    Relaxed.Standard.MuH.CD8$AICc,  Relaxed.Standard.MuH.CID2$AICc,   
                                    Relaxed.Standard.MuH.CID3$AICc, Relaxed.Standard.MuH.CID4$AICc,   
                                    Relaxed.Standard.MuH.CID5$AICc ,Relaxed.Standard.MuH.CID6$AICc,   
                                    Relaxed.Standard.MuH.CID7$AICc, Relaxed.Standard.MuH.CID8$AICc,   
                                    Relaxed.Standard.Null.MuS$AICc, Relaxed.Standard.True.MuS$AICc)

LogL_Relaxed.Standard_SSE_Full <- c(Relaxed.Standard.MuH.CD2$loglik,  Relaxed.Standard.MuH.CD3$loglik ,   
                                    Relaxed.Standard.MuH.CD4$loglik,  Relaxed.Standard.MuH.CD5$loglik,    
                                    Relaxed.Standard.MuH.CD6$loglik,  Relaxed.Standard.MuH.CD7$loglik,    
                                    Relaxed.Standard.MuH.CD8$loglik,  Relaxed.Standard.MuH.CID2$loglik,   
                                    Relaxed.Standard.MuH.CID3$loglik, Relaxed.Standard.MuH.CID4$loglik,   
                                    Relaxed.Standard.MuH.CID5$loglik ,Relaxed.Standard.MuH.CID6$loglik,   
                                    Relaxed.Standard.MuH.CID7$loglik, Relaxed.Standard.MuH.CID8$loglik,   
                                    Relaxed.Standard.Null.MuS$loglik, Relaxed.Standard.True.MuS$loglik)


Res_Relaxed.Standard_AICc_SSE_Full <-as.data.frame(AICc_Relaxed.Standard_SSE_Full, row.names = c("Relaxed.Standard.MuH.CD2",  "Relaxed.Standard.MuH.CD3" ,   
                                                                                                 "Relaxed.Standard.MuH.CD4",  "Relaxed.Standard.MuH.CD5",    
                                                                                                 "Relaxed.Standard.MuH.CD6",  "Relaxed.Standard.MuH.CD7",    
                                                                                                 "Relaxed.Standard.MuH.CD8",  "Relaxed.Standard.MuH.CID2",   
                                                                                                 "Relaxed.Standard.MuH.CID3", "Relaxed.Standard.MuH.CID4",   
                                                                                                 "Relaxed.Standard.MuH.CID5" ,"Relaxed.Standard.MuH.CID6",   
                                                                                                 "Relaxed.Standard.MuH.CID7", "Relaxed.Standard.MuH.CID8",   
                                                                                                 "Relaxed.Standard.Null.MuS", "Relaxed.Standard.True.MuS"))

Res_Relaxed.Standard_LogL_SSE_Full <-as.data.frame(LogL_Relaxed.Standard_SSE_Full, row.names = c("Relaxed.Standard.MuH.CD2",  "Relaxed.Standard.MuH.CD3" ,   
                                                                                                 "Relaxed.Standard.MuH.CD4",  "Relaxed.Standard.MuH.CD5",    
                                                                                                 "Relaxed.Standard.MuH.CD6",  "Relaxed.Standard.MuH.CD7",    
                                                                                                 "Relaxed.Standard.MuH.CD8",  "Relaxed.Standard.MuH.CID2",   
                                                                                                 "Relaxed.Standard.MuH.CID3", "Relaxed.Standard.MuH.CID4",   
                                                                                                 "Relaxed.Standard.MuH.CID5" ,"Relaxed.Standard.MuH.CID6",   
                                                                                                 "Relaxed.Standard.MuH.CID7", "Relaxed.Standard.MuH.CID8",   
                                                                                                 "Relaxed.Standard.Null.MuS", "Relaxed.Standard.True.MuS"))

Res_AIC_Relaxed.Standard_SSE_Full <-as.data.frame(AIC_Relaxed.Standard_SSE_Full, row.names = c("Relaxed.Standard.MuH.CD2",  "Relaxed.Standard.MuH.CD3" ,   
                                                                                               "Relaxed.Standard.MuH.CD4",  "Relaxed.Standard.MuH.CD5",    
                                                                                               "Relaxed.Standard.MuH.CD6",  "Relaxed.Standard.MuH.CD7",    
                                                                                               "Relaxed.Standard.MuH.CD8",  "Relaxed.Standard.MuH.CID2",   
                                                                                               "Relaxed.Standard.MuH.CID3", "Relaxed.Standard.MuH.CID4",   
                                                                                               "Relaxed.Standard.MuH.CID5" ,"Relaxed.Standard.MuH.CID6",   
                                                                                               "Relaxed.Standard.MuH.CID7", "Relaxed.Standard.MuH.CID8",   
                                                                                               "Relaxed.Standard.Null.MuS", "Relaxed.Standard.True.MuS"))

Res_Relaxed.Standard_LogL_SSE_Full$AICc<-Res_Relaxed.Standard_AICc_SSE_Full$AICc_Relaxed.Standard_SSE_Full #Strict_SYM_NoDual_HiddenRatCat2_Model wins at 98.2% wAIC
Res_Relaxed.Standard_LogL_SSE_Full$AIC <- Res_AIC_Relaxed.Standard_SSE_Full$AIC_Relaxed.Standard_SSE_Full
Res_Relaxed.Standard_LogL_SSE_Full$AICw <- aicw(AICc_Relaxed.Standard_SSE_Full)
#Write.csv if you wish

## Transition rates

RelaxedStandardSolution <- Relaxed.Standard.MuH.CID8$solution

## Transpose the dataframe
RelaxedStandardSolution.Transposed <- as.data.frame(t(as.matrix(RelaxedStandardSolution)))

## a little help from dplyr to keep only the columns with underscores (transition rate estimates)
RelaxedStandardSolution.Rates <- RelaxedStandardSolution.Transposed %>% select(contains("_"))

## 1
RelaxedStandardSolution.Rates1 <- RelaxedStandardSolution.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
RelaxedStandardSolution.Rates.00.10 <- RelaxedStandardSolution.Rates1 %>% select(contains("q00"))
## q10 > q00
RelaxedStandardSolution.Rates.10.00 <- RelaxedStandardSolution.Rates1 %>% select(contains("q10"))

## 2
RelaxedStandardSolution.Rates.q01q00 <- RelaxedStandardSolution.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
RelaxedStandardSolution.Rates.00.01 <- RelaxedStandardSolution.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
RelaxedStandardSolution.Rates.01.00 <- RelaxedStandardSolution.Rates.q01q00 %>% select(contains("q01"))

## 3
RelaxedStandardSolution.Rates.q01q11 <- RelaxedStandardSolution.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
RelaxedStandardSolution.Rates.01.11 <- RelaxedStandardSolution.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
RelaxedStandardSolution.Rates.11.01 <- RelaxedStandardSolution.Rates.q01q11 %>% select(contains("q11"))

## 4
RelaxedStandardSolution.Rates.q10q11 <- RelaxedStandardSolution.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
RelaxedStandardSolution.Rates.10.11 <- RelaxedStandardSolution.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
RelaxedStandardSolution.Rates.11.10 <- RelaxedStandardSolution.Rates.q10q11 %>% select(contains("q11"))

## Combine all into a new dataframe for plotting
RelaxedStandard_qrates <- list(RelaxedStandardSolution.Rates.00.10, 
                               RelaxedStandardSolution.Rates.10.00, 
                               RelaxedStandardSolution.Rates.00.01, 
                               RelaxedStandardSolution.Rates.01.00, 
                               RelaxedStandardSolution.Rates.01.11, 
                               RelaxedStandardSolution.Rates.11.01, 
                               RelaxedStandardSolution.Rates.10.11, 
                               RelaxedStandardSolution.Rates.11.10) %>%
  setNames(c("00.10",
             "10.00",
             "00.01",
             "01.00",
             "01.11",
             "11.01",
             "10.11",
             "11.10"))



#Just sample one
for(i in 1:length(RelaxedStandard_qrates)){
  MuHiSSERelaxedStandardRates <- sample(RelaxedStandard_qrates[[i]][1]) * 100
  MuHiSSERelaxedStandardRatesPrint <- print(round(MuHiSSERelaxedStandardRates, 4))
  #write.csv(MuHiSSERelaxedStandardRates[[i]][1], "MuHiSSERelaxedStandardRates.Strict.MuH.CID8.csv")
}

## Reconstruct ancestral states for this MuHCID8

Relaxed.Standard.MuH.CID8.recon <-
  hisse::MarginReconMuHiSSE(
    phy = Relaxed.Standard.MuH.CID8$phy,
    data = Relaxed.Standard.MuH.CID8$data,
    f = Relaxed.Standard.MuH.CID8$f,
    pars = Relaxed.Standard.MuH.CID8$solution,
    hidden.states = 8,
    aic = Relaxed.Standard.MuH.CID8$AICc,
    n.cores = 3,
    verbose = TRUE
  )

## Plotting
Relaxed.Standard.MuH.CID8.recon 
Relaxed.Standard.MuH.CID8.recon.proc <- m_process_recon(Relaxed.Standard.MuH.CID8.recon)
Relaxed.Standard.MuH.CID8.recon.recon.AncestralStates <- m_trait_recon(
  processed_recon = Relaxed.Standard.MuH.CID8.recon.proc,
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
Relaxed.Standard.MuH.CID8.recon.recon.AncestralStates + geom_cladelabel(node=2362, label="Gekkota", align=T, geom='label', fill='white',offset.text = 6)+
  geom_cladelabel(node=3052, label="Scincidae", align=T, geom='label', fill='white',offset.text = 40)+
  geom_cladelabel(node=4526, label="Dactyloidae", align=T, geom='label', fill='white',offset.text = 6)

## Net Diverisification
RelaxedStandard.MuH.CID8.NetDiv <- m_scatterplot(
  processed_recon = Relaxed.Standard.MuH.CID8.recon.proc,
  parameter = 'net.div',
  states_names = c('00','01','10','11'), colors = c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"), 
  plot_as_waiting_time = FALSE) + theme_classic()

RelaxedStandard.MuH.CID8.NetDiv + ylab("Net Diversfication")
dev.size("in") #7.375000 5.305556

