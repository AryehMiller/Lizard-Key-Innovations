## fitPagel with fixed ("fossilized") root state
## Full and Standard (referred to here as "limbless") datasets

library(phytools)
library(ape)
library(data.table)
library(ggplot2)
library(viridis)
library(Hmisc)
library(ggpubr)

##################################
######                      ######
######   Strict "Standard"  ######
######                      ######
##################################

## Read in trait dataset

LimblessDataset<- read.csv("Meiri2360_StandardEcologyLimbed.csv", row.names = 1)

## Read in phylogenies
LimblessTree <- read.nexus(file = "Meiri23Dec20_LimblessRemoved.tre")

  
## Add tip of zero-length to the root using "bind.tip" and then assigning it non-arboreal and padless
Limbless.tree.root <-bind.tip(LimblessTree,tip.label="ROOT",edge.length=0,
                          where=Ntip(LimblessTree)+1)
#plotTree(Limbless.tree.root,ftype="off",lwd=1)
nnL<-which(Limbless.tree.root$tip.label=="ROOT")
tiplabels(Limbless.tree.root$tip.label[nnL],nnL,adj=c(-0.1,0.5),
          frame="none",cex=0.8,font=3)

## Add fossilized root states
arboreal.strict.Limbless<-as.factor(setNames(LimblessDataset[,4],rownames(LimblessDataset)))
toepads.Limbless<-as.factor(setNames(LimblessDataset[,7],rownames(LimblessDataset)))

Limbless.arboreal.strict.fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                  setNames(as.character(arboreal.strict.Limbless),names(arboreal.strict.Limbless))))
Limbless.padless.fossilized <-as.factor(c(setNames("0","ROOT"),
                                 setNames(as.character(toepads.Limbless),names(toepads.Limbless))))

#fit fixed-root models for Strict Limbless dataset

##ARD

## Dependent y model
strict.Limbless.fit.fixed.dep.y.ARD<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="y",model="ARD") #strict

## Dependent x model
strict.Limbless.fit.fixed.dep.x.ARD <-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="x",model="ARD") #strict

## Interdependent model
strict.Limbless.fit.fixed.dep.xy.ARD<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="xy",model="ARD") #strict

## SYM

## Dependent y model
strict.Limbless.fit.fixed.dep.y.SYM<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="y",model="SYM") #strict

## Dependent x model
strict.Limbless.fit.fixed.dep.x.SYM <-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="x",model="SYM") #strict

## Interdependent model
strict.Limbless.fit.fixed.dep.xy.SYM<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="xy",model="SYM") #strict

## ER

## Dependent y model
strict.Limbless.fit.fixed.dep.y.ER<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="y",model="ER") #strict

## Dependent x model
strict.Limbless.fit.fixed.dep.x.ER <-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="x",model="ER") #strict

## Interdependent model
strict.Limbless.fit.fixed.dep.xy.ER<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.strict.fossilized,dep.var="xy",model="ER") #strict

strict.Limbless.aic <-setNames(c(strict.Limbless.fit.fixed.dep.y.ARD$dependent.AIC,
                                  strict.Limbless.fit.fixed.dep.x.ARD$dependent.AIC,
                                  strict.Limbless.fit.fixed.dep.xy.ARD$dependent.AIC,
                                  strict.Limbless.fit.fixed.dep.y.ER$dependent.AIC,
                                  strict.Limbless.fit.fixed.dep.x.ER$dependent.AIC,
                                  strict.Limbless.fit.fixed.dep.xy.ER$dependent.AIC),
                                c("ARD dependent y",
                                  "ARD dependent x",
                                  "ARD dependent x&y",
                                  "ER dependent y",
                                  "ER dependent x",
                                  "ER dependent x&y"))

print(strict.Limbless.aic)
aicw(strict.Limbless.aic)

plot(strict.Limbless.fit.fixed.dep.xy.ARD)
##################################
######                      ######
######  Relaxed "Limbless"  ######
######                      ######
##################################

## Add fossilized root states
arboreal.relaxed.Limbless<-as.factor(setNames(LimblessDataset[,5],rownames(LimblessDataset)))

Limbless.arboreal.relaxed.fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                                  setNames(as.character(arboreal.relaxed.Limbless),names(arboreal.relaxed.Limbless))))

## Fit fixed-root models for Relaxed Limbless dataset

## ARD

## Dependent y model
relaxed.Limbless.fit.fixed.dep.y.ARD<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="y",model="ARD") #strict

#relaxed.Limbless.fit.fixed.dep.y.ARD
#plot(relaxed.Limbless.fit.fixed.dep.y.ARD)
#title("relaxed.Limbless.fit.fixed.dep.y.ARD")

## Dependent x model
relaxed.Limbless.fit.fixed.dep.x.ARD <-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="x",model="ARD") #strict

## Interdependent model
relaxed.Limbless.fit.fixed.dep.xy.ARD<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="xy",model="ARD") #strict

## SYM

## Dependent y model
relaxed.Limbless.fit.fixed.dep.y.SYM<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="y",model="SYM") #strict

## Dependent x model
relaxed.Limbless.fit.fixed.dep.x.SYM <-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="x",model="SYM") #strict

## Interdependent model
relaxed.Limbless.fit.fixed.dep.xy.SYM<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="xy",model="SYM") #strict

## ER

## Dependent y model
relaxed.Limbless.fit.fixed.dep.y.ER<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="y",model="ER") #strict

## Dependent x model
relaxed.Limbless.fit.fixed.dep.x.ER <-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="x",model="ER") #strict

## Interdependent model
relaxed.Limbless.fit.fixed.dep.xy.ER<-fitPagel(Limbless.tree.root,Limbless.padless.fossilized,Limbless.arboreal.relaxed.fossilized,dep.var="xy",model="ER") #strict

## Organize for AIC scores

relaxed.Limbless.ER.aic <-setNames(c(relaxed.Limbless.fit.fixed.dep.y.ER$dependent.AIC,
                                    relaxed.Limbless.fit.fixed.dep.x.ER$dependent.AIC,
                                    relaxed.Limbless.fit.fixed.dep.xy.ER$dependent.AIC),
                                  c("dependent y",
                                    "dependent x","dependent x&y"))
relaxed.Limbless.ER.aic
aic.w(relaxed.Limbless.ER.aic)


relaxed.Limbless.aic <-setNames(c(relaxed.Limbless.fit.fixed.dep.y.ARD$dependent.AIC,
                                     relaxed.Limbless.fit.fixed.dep.x.ARD$dependent.AIC,
                                     relaxed.Limbless.fit.fixed.dep.xy.ARD$dependent.AIC,
                                     relaxed.Limbless.fit.fixed.dep.y.ER$dependent.AIC,
                                     relaxed.Limbless.fit.fixed.dep.x.ER$dependent.AIC,
                                     relaxed.Limbless.fit.fixed.dep.xy.ER$dependent.AIC),
                                   c("ARD dependent y",
                                     "ARD dependent x",
                                     "ARD dependent x&y",
                                     "ER dependent y",
                                     "ER dependent x",
                                     "ER dependent x&y"))

print(relaxed.Limbless.aic)
aicw(relaxed.Limbless.aic)

plot(relaxed.Limbless.fit.fixed.dep.y.ARD)
##################################
######                      ######
######   Strict "Full"      ######
######                      ######
##################################

#Full trait dataset
FullDataset<- read.csv("Meiri_Pyron_2692_5Dec20_ForAnalysis.csv", row.names = 1)

#Full phylogeny
FullTree <- read.nexus(file = "pruned.tree.Meiri.RevisedPyron.tre")

## Add tip of zero-length to the root using "bind.tip" and then assigning it non-arboreal and padless
Full.tree.root <-bind.tip(FullTree,tip.label="ROOT",edge.length=0,
                              where=Ntip(FullTree)+1)
#plotTree(Full.tree.root,ftype="off",lwd=1)
nnF<-which(Full.tree.root$tip.label=="ROOT")
tiplabels(Full.tree.root$tip.label[nnF],nnF,adj=c(-0.1,0.5),
          frame="none",cex=0.8,font=3)

## Add fossilized root states
arboreal.strict.Full<-as.factor(setNames(FullDataset[,4],rownames(FullDataset)))
toepads.Full<-as.factor(setNames(FullDataset[,7],rownames(FullDataset)))

Full.arboreal.strict.fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                                  setNames(as.character(arboreal.strict.Full),names(arboreal.strict.Full))))
Full.padless.fossilized <-as.factor(c(setNames("0","ROOT"),
                                          setNames(as.character(toepads.Full),names(toepads.Full))))

#fit fixed-root models for Strict Limbless dataset

##ARD

## Dependent y model
strict.Full.fit.fixed.dep.y.ARD<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="y",model="ARD") #strict

## Dependent x model
strict.Full.fit.fixed.dep.x.ARD <-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="x",model="ARD") #strict

## Interdependent model
strict.Full.fit.fixed.dep.xy.ARD<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="xy",model="ARD") #strict

## SYM

## Dependent y model
strict.Full.fit.fixed.dep.y.SYM<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="y",model="SYM") #strict

## Dependent x model
strict.Full.fit.fixed.dep.x.SYM <-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="x",model="SYM") #strict

## Interdependent model
strict.Full.fit.fixed.dep.xy.SYM<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="xy",model="SYM") #strict

## ER

## Dependent y model
strict.Full.fit.fixed.dep.y.ER<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="y",model="ER") #strict

## Dependent x model
strict.Full.fit.fixed.dep.x.ER <-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="x",model="ER") #strict

## Interdependent model
strict.Full.fit.fixed.dep.xy.ER<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="xy",model="ER") #strict

strict.Full.aic <-setNames(c(strict.Full.fit.fixed.dep.y.ARD$dependent.AIC,
                                 strict.Full.fit.fixed.dep.x.ARD$dependent.AIC,
                                 strict.Full.fit.fixed.dep.xy.ARD$dependent.AIC,
                                 strict.Full.fit.fixed.dep.y.ER$dependent.AIC,
                                 strict.Full.fit.fixed.dep.x.ER$dependent.AIC,
                                 strict.Full.fit.fixed.dep.xy.ER$dependent.AIC),
                               c("ARD dependent y",
                                 "ARD dependent x",
                                 "ARD dependent x&y",
                                 "ER dependent y",
                                 "ER dependent x",
                                 "ER dependent x&y"))

print(strict.Full.aic)
aicw(strict.Full.aic)
plot(strict.Full.fit.fixed.dep.xy.ARD) #Best Strict Full model
##################################
######                      ######
######   Relaxed "Full"     ######
######                      ######
##################################

## Add fossilized root states
arboreal.relaxed.Full<-as.factor(setNames(FullDataset[,5],rownames(FullDataset)))

Full.arboreal.relaxed.fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                              setNames(as.character(arboreal.relaxed.Full),names(arboreal.relaxed.Full))))

#fit fixed-root models for Relaxed Full dataset

##ARD

## Dependent y model
relaxed.Full.fit.fixed.dep.y.ARD<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="y",model="ARD") #relaxed

## Dependent x model
relaxed.Full.fit.fixed.dep.x.ARD <-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="x",model="ARD") #relaxed

## Interdependent model
relaxed.Full.fit.fixed.dep.xy.ARD<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="xy",model="ARD") #relaxed

## SYM

## Dependent y model
relaxed.Full.fit.fixed.dep.y.SYM<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="y",model="SYM") #relaxed

## Dependent x model
relaxed.Full.fit.fixed.dep.x.SYM <-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="x",model="SYM") #relaxed

## Interdependent model
relaxed.Full.fit.fixed.dep.xy.SYM<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="xy",model="SYM") #relaxed

## ER

## Dependent y model
relaxed.Full.fit.fixed.dep.y.ER<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="y",model="ER") #relaxed

## Dependent x model
relaxed.Full.fit.fixed.dep.x.ER <-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="x",model="ER") #relaxed

## Interdependent model
relaxed.Full.fit.fixed.dep.xy.ER<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.relaxed.fossilized,dep.var="xy",model="ER") #relaxed

relaxed.Full.aic <-setNames(c(relaxed.Full.fit.fixed.dep.y.ARD$dependent.AIC,
                             relaxed.Full.fit.fixed.dep.x.ARD$dependent.AIC,
                             relaxed.Full.fit.fixed.dep.xy.ARD$dependent.AIC,
                             relaxed.Full.fit.fixed.dep.y.ER$dependent.AIC,
                             relaxed.Full.fit.fixed.dep.x.ER$dependent.AIC,
                             relaxed.Full.fit.fixed.dep.xy.ER$dependent.AIC),
                           c("ARD dependent y",
                             "ARD dependent x",
                             "ARD dependent x&y",
                             "ER dependent y",
                             "ER dependent x",
                             "ER dependent x&y"))

print(relaxed.Full.aic)
aicw(relaxed.Full.aic)
plot(relaxed.Full.fit.fixed.dep.xy.ARD) #Best Relaxed Full Model

## Plot results with ggplot2
## Gather rates

#Relaxed Full

relaxed.Full.fit.fixed.dep.y.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Full.fit.fixed.dep.y.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Full.fit.fixed.dep.x.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Full.fit.fixed.dep.x.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Full.fit.fixed.dep.xy.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Full.fit.fixed.dep.xy.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Full.fit.fixed.dep.y.ER.Q <- as.data.frame(t(as.matrix(relaxed.Full.fit.fixed.dep.y.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Full.fit.fixed.dep.x.ER.Q <- as.data.frame(t(as.matrix(relaxed.Full.fit.fixed.dep.x.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Full.fit.fixed.dep.xy.ER.Q <- as.data.frame(t(as.matrix(relaxed.Full.fit.fixed.dep.xy.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))

##To export table with labels
#relaxed.Full.fit.fixed.dep.y.ARD.Q <- setnames(relaxed.Full.fit.fixed.dep.y.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Full.fit.fixed.dep.x.ARD.Q <- setnames(relaxed.Full.fit.fixed.dep.x.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Full.fit.fixed.dep.xy.ARD.Q <- setnames(relaxed.Full.fit.fixed.dep.xy.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Full.fit.fixed.dep.y.ER.Q <- setnames(relaxed.Full.fit.fixed.dep.y.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Full.fit.fixed.dep.x.ER.Q <- setnames(relaxed.Full.fit.fixed.dep.x.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Full.fit.fixed.dep.xy.ER.Q <- setnames(relaxed.Full.fit.fixed.dep.xy.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#
#Relaxed.Full.Q.Export <-as.data.frame(Relaxed.Full.Q, row.names = c("ARD dependent y",
#                                                             "ARD dependent x",
#                                                             "ARD dependent x&y",
#                                                             "ER dependent y",
#                                                             "ER dependent x",
#                                                             "ER dependent x&y"))
#
#Relaxed.Full.AIC <- c(relaxed.Full.fit.fixed.dep.y.ARD$dependent.AIC,
#                      relaxed.Full.fit.fixed.dep.x.ARD$dependent.AIC,
#                      relaxed.Full.fit.fixed.dep.xy.ARD$dependent.AIC,
#                      relaxed.Full.fit.fixed.dep.y.ER$dependent.AIC,
#                      relaxed.Full.fit.fixed.dep.x.ER$dependent.AIC,
#                      relaxed.Full.fit.fixed.dep.xy.ER$dependent.AIC)
#
#
#Relaxed.Full.Q.Export$AIC<-Relaxed.Full.AIC 

###

relaxed.Full.fit.fixed.dep.y.ARD.Q <- setnames(relaxed.Full.fit.fixed.dep.y.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Full.fit.fixed.dep.x.ARD.Q <- setnames(relaxed.Full.fit.fixed.dep.x.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Full.fit.fixed.dep.xy.ARD.Q <- setnames(relaxed.Full.fit.fixed.dep.xy.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Full.fit.fixed.dep.y.ER.Q <- setnames(relaxed.Full.fit.fixed.dep.y.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Full.fit.fixed.dep.x.ER.Q <- setnames(relaxed.Full.fit.fixed.dep.x.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Full.fit.fixed.dep.xy.ER.Q <- setnames(relaxed.Full.fit.fixed.dep.xy.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))

Relaxed.Full.Q <- rbind(relaxed.Full.fit.fixed.dep.y.ARD.Q,
                relaxed.Full.fit.fixed.dep.x.ARD.Q,
                relaxed.Full.fit.fixed.dep.xy.ARD.Q,
                relaxed.Full.fit.fixed.dep.y.ER.Q,
                relaxed.Full.fit.fixed.dep.x.ER.Q,
                relaxed.Full.fit.fixed.dep.xy.ER.Q)

Relaxed.Full.Q <-as.data.frame(Relaxed.Full.Q, row.names = c("ARD dependent y",
                                                             "ARD dependent x",
                                                             "ARD dependent x&y",
                                                             "ER dependent y",
                                                             "ER dependent x",
                                                             "ER dependent x&y"))


## Add in AIC scores
Relaxed.Full.Q.Export$AIC<-Relaxed.Full.AIC 

## Finalize table
setDT(Relaxed.Full.Q, keep.rownames = TRUE)
colnames(Relaxed.Full.Q)[1] <- "Model"

library(reshape2)
Relaxed.Full.Q.Melted <- melt((Relaxed.Full.Q))

Relaxed.Full.Q.Melted.level.order <- factor(Relaxed.Full.Q.Melted$variable, level = c('PadlessNonArboreal.PadbearingNonArboreal',
                                                                                    'PadlessArboreal.PadbearingArboreal',
                                                                                    'PadbearingArboreal.PadlessArboreal',
                                                                                    'PadbearingNonArboreal.PadlessNonArboreal',
                                                                                    'PadbearingArboreal.PadbearingNonArboreal',
                                                                                    'PadlessNonArboreal.PadlessArboreal',
                                                                                    'PadlessArboreal.PadlessNonArboreal', 
                                                                                    'PadbearingNonArboreal.PadbearingArboreal'))


Relaxed.Full.Q.Jitter <- ggplot(Relaxed.Full.Q.Melted, aes(x = Relaxed.Full.Q.Melted.level.order, y = value))+ geom_jitter(aes(color = Model), 
                position = position_jitter(0), size = 2.5) +
  xlab("") + ylab("Transition Rate")+ggtitle("Relaxed Full Dataset: Pagel")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_colour_viridis_d(option = "viridis")+theme_classic()+coord_flip(ylim = c(0,0.03))

dev.size("in") #[1] 10.77778  5.87500

#Strict Full
strict.Full.fit.fixed.dep.y.ARD.Q <- as.data.frame(t(as.matrix(strict.Full.fit.fixed.dep.y.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Full.fit.fixed.dep.x.ARD.Q <- as.data.frame(t(as.matrix(strict.Full.fit.fixed.dep.x.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Full.fit.fixed.dep.xy.ARD.Q <- as.data.frame(t(as.matrix(strict.Full.fit.fixed.dep.xy.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Full.fit.fixed.dep.y.ER.Q <- as.data.frame(t(as.matrix(strict.Full.fit.fixed.dep.y.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Full.fit.fixed.dep.x.ER.Q <- as.data.frame(t(as.matrix(strict.Full.fit.fixed.dep.x.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Full.fit.fixed.dep.xy.ER.Q <- as.data.frame(t(as.matrix(strict.Full.fit.fixed.dep.xy.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))

##Export table with labels
#strict.Full.fit.fixed.dep.y.ARD.Q <- setnames(strict.Full.fit.fixed.dep.y.ARD.Q,  c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Full.fit.fixed.dep.x.ARD.Q <- setnames(strict.Full.fit.fixed.dep.x.ARD.Q,  c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Full.fit.fixed.dep.xy.ARD.Q <- setnames(strict.Full.fit.fixed.dep.xy.ARD.Q,  c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Full.fit.fixed.dep.y.ER.Q <- setnames(strict.Full.fit.fixed.dep.y.ER.Q,  c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Full.fit.fixed.dep.x.ER.Q <- setnames(strict.Full.fit.fixed.dep.x.ER.Q,  c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Full.fit.fixed.dep.xy.ER.Q <- setnames(strict.Full.fit.fixed.dep.xy.ER.Q,  c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#
#strict.Full.Q.export <- rbind(strict.Full.fit.fixed.dep.y.ARD.Q,
#                       strict.Full.fit.fixed.dep.x.ARD.Q,
#                       strict.Full.fit.fixed.dep.xy.ARD.Q,
#                       strict.Full.fit.fixed.dep.y.ER.Q,
#                       strict.Full.fit.fixed.dep.x.ER.Q,
#                       strict.Full.fit.fixed.dep.xy.ER.Q)
#
#strict.Full.Q.export <-as.data.frame(strict.Full.Q.export, row.names = c("ARD dependent y",
#                                                           "ARD dependent x",
#                                                           "ARD dependent x&y",
#                                                           "ER dependent y",
#                                                           "ER dependent x",
#                                                           "ER dependent x&y"))
#
#strict.Full.AIC <- c(strict.Full.fit.fixed.dep.y.ARD$dependent.AIC,
#                     strict.Full.fit.fixed.dep.x.ARD$dependent.AIC,
#                     strict.Full.fit.fixed.dep.xy.ARD$dependent.AIC,
#                     strict.Full.fit.fixed.dep.y.ER$dependent.AIC,
#                     strict.Full.fit.fixed.dep.x.ER$dependent.AIC,
#                     strict.Full.fit.fixed.dep.xy.ER$dependent.AIC)
#
### Add in AIC scores
#strict.Full.Q.export$AIC<-strict.Full.AIC 
#
#strict.Full.Q.export

####
strict.Full.fit.fixed.dep.y.ARD.Q <- setnames(strict.Full.fit.fixed.dep.y.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Full.fit.fixed.dep.x.ARD.Q <- setnames(strict.Full.fit.fixed.dep.x.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Full.fit.fixed.dep.xy.ARD.Q <- setnames(strict.Full.fit.fixed.dep.xy.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Full.fit.fixed.dep.y.ER.Q <- setnames(strict.Full.fit.fixed.dep.y.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Full.fit.fixed.dep.x.ER.Q <- setnames(strict.Full.fit.fixed.dep.x.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Full.fit.fixed.dep.xy.ER.Q <- setnames(strict.Full.fit.fixed.dep.xy.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))

strict.Full.Q <- rbind(strict.Full.fit.fixed.dep.y.ARD.Q,
                        strict.Full.fit.fixed.dep.x.ARD.Q,
                        strict.Full.fit.fixed.dep.xy.ARD.Q,
                        strict.Full.fit.fixed.dep.y.ER.Q,
                        strict.Full.fit.fixed.dep.x.ER.Q,
                        strict.Full.fit.fixed.dep.xy.ER.Q)

strict.Full.Q <-as.data.frame(strict.Full.Q, row.names = c("ARD dependent y",
                                                             "ARD dependent x",
                                                             "ARD dependent x&y",
                                                             "ER dependent y",
                                                             "ER dependent x",
                                                             "ER dependent x&y"))


## Finalize table
setDT(strict.Full.Q, keep.rownames = TRUE)
colnames(strict.Full.Q)[1] <- "Model"

library(reshape2)
strict.Full.Q.Melted <- melt((strict.Full.Q))

strict.Full.Q.Melted.level.order <- factor(strict.Full.Q.Melted$variable, level = c('PadlessNonArboreal.PadbearingNonArboreal',
                                                                                            'PadlessArboreal.PadbearingArboreal',
                                                                                            'PadbearingArboreal.PadlessArboreal',
                                                                                            'PadbearingNonArboreal.PadlessNonArboreal',
                                                                                            'PadbearingArboreal.PadbearingNonArboreal',
                                                                                            'PadlessNonArboreal.PadlessArboreal',
                                                                                            'PadlessArboreal.PadlessNonArboreal', 
                                                                                            'PadbearingNonArboreal.PadbearingArboreal'))


Strict.Full.Q.Jitter <- ggplot(strict.Full.Q.Melted, aes(x = strict.Full.Q.Melted.level.order, y = value))+ geom_jitter(aes(color = Model), 
                                                                         position = position_jitter(0), size = 2.5) +
  xlab("") + ylab("Transition Rate")+ggtitle("Strict Full Dataset: Pagel")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_colour_viridis_d(option = "viridis")+theme_classic()+coord_flip(ylim = c(0,0.03))

dev.size("in") #[1] 10.77778  5.87500

#Strict Standard

strict.Limbless.fit.fixed.dep.y.ARD.Q <- as.data.frame(t(as.matrix(strict.Limbless.fit.fixed.dep.y.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Limbless.fit.fixed.dep.x.ARD.Q <- as.data.frame(t(as.matrix(strict.Limbless.fit.fixed.dep.x.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Limbless.fit.fixed.dep.xy.ARD.Q <- as.data.frame(t(as.matrix(strict.Limbless.fit.fixed.dep.xy.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Limbless.fit.fixed.dep.y.ER.Q <- as.data.frame(t(as.matrix(strict.Limbless.fit.fixed.dep.y.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Limbless.fit.fixed.dep.x.ER.Q <- as.data.frame(t(as.matrix(strict.Limbless.fit.fixed.dep.x.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Limbless.fit.fixed.dep.xy.ER.Q <- as.data.frame(t(as.matrix(strict.Limbless.fit.fixed.dep.xy.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))

### Table for export

#strict.Limbless.fit.fixed.dep.y.ARD.Q <- setnames(strict.Limbless.fit.fixed.dep.y.ARD.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Limbless.fit.fixed.dep.x.ARD.Q <- setnames(strict.Limbless.fit.fixed.dep.x.ARD.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Limbless.fit.fixed.dep.xy.ARD.Q <- setnames(strict.Limbless.fit.fixed.dep.xy.ARD.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Limbless.fit.fixed.dep.y.ER.Q <- setnames(strict.Limbless.fit.fixed.dep.y.ER.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Limbless.fit.fixed.dep.x.ER.Q <- setnames(strict.Limbless.fit.fixed.dep.x.ER.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Limbless.fit.fixed.dep.xy.ER.Q <- setnames(strict.Limbless.fit.fixed.dep.xy.ER.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#
#strict.Limbless.Q.export <- rbind(strict.Limbless.fit.fixed.dep.y.ARD.Q,
#                           strict.Limbless.fit.fixed.dep.x.ARD.Q,
#                           strict.Limbless.fit.fixed.dep.xy.ARD.Q,
#                           strict.Limbless.fit.fixed.dep.y.ER.Q,
#                           strict.Limbless.fit.fixed.dep.x.ER.Q,
#                           strict.Limbless.fit.fixed.dep.xy.ER.Q)
#
#strict.Limbless.Q.export <-as.data.frame(strict.Limbless.Q.export, row.names = c("ARD dependent y",
#                                                                   "ARD dependent x",
#                                                                   "ARD dependent x&y",
#                                                                   "ER dependent y",
#                                                                   "ER dependent x",
#                                                                   "ER dependent x&y"))
#
#strict.Limbless.AIC <- c(strict.Limbless.fit.fixed.dep.y.ARD$dependent.AIC,
#                         strict.Limbless.fit.fixed.dep.x.ARD$dependent.AIC,
#                         strict.Limbless.fit.fixed.dep.xy.ARD$dependent.AIC,
#                         strict.Limbless.fit.fixed.dep.y.ER$dependent.AIC,
#                         strict.Limbless.fit.fixed.dep.x.ER$dependent.AIC,
#                         strict.Limbless.fit.fixed.dep.xy.ER$dependent.AIC)
#
### Add in AIC scores
#strict.Limbless.Q.export$AIC<-strict.Limbless.AIC 
#
#strict.Limbless.Q.export

###

strict.Limbless.fit.fixed.dep.y.ARD.Q <- setnames(strict.Limbless.fit.fixed.dep.y.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Limbless.fit.fixed.dep.x.ARD.Q <- setnames(strict.Limbless.fit.fixed.dep.x.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Limbless.fit.fixed.dep.xy.ARD.Q <- setnames(strict.Limbless.fit.fixed.dep.xy.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Limbless.fit.fixed.dep.y.ER.Q <- setnames(strict.Limbless.fit.fixed.dep.y.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Limbless.fit.fixed.dep.x.ER.Q <- setnames(strict.Limbless.fit.fixed.dep.x.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Limbless.fit.fixed.dep.xy.ER.Q <- setnames(strict.Limbless.fit.fixed.dep.xy.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))

strict.Limbless.Q <- rbind(strict.Limbless.fit.fixed.dep.y.ARD.Q,
                       strict.Limbless.fit.fixed.dep.x.ARD.Q,
                       strict.Limbless.fit.fixed.dep.xy.ARD.Q,
                       strict.Limbless.fit.fixed.dep.y.ER.Q,
                       strict.Limbless.fit.fixed.dep.x.ER.Q,
                       strict.Limbless.fit.fixed.dep.xy.ER.Q)

strict.Limbless.Q <-as.data.frame(strict.Limbless.Q, row.names = c("ARD dependent y",
                                                           "ARD dependent x",
                                                           "ARD dependent x&y",
                                                           "ER dependent y",
                                                           "ER dependent x",
                                                           "ER dependent x&y"))



## Finalize table
setDT(strict.Limbless.Q, keep.rownames = TRUE)
colnames(strict.Limbless.Q)[1] <- "Model"

library(reshape2)
strict.Limbless.Q.Melted <- melt((strict.Limbless.Q))

strict.Limbless.Q.Melted.level.order <- factor(strict.Limbless.Q.Melted$variable, level = c('PadlessNonArboreal.PadbearingNonArboreal',
                                                                                              'PadlessArboreal.PadbearingArboreal',
                                                                                              'PadbearingArboreal.PadlessArboreal',
                                                                                              'PadbearingNonArboreal.PadlessNonArboreal',
                                                                                              'PadbearingArboreal.PadbearingNonArboreal',
                                                                                              'PadlessNonArboreal.PadlessArboreal',
                                                                                              'PadlessArboreal.PadlessNonArboreal', 
                                                                                              'PadbearingNonArboreal.PadbearingArboreal'))

Strict.Standard.Q.Jitter <-ggplot(strict.Limbless.Q.Melted, aes(x = strict.Limbless.Q.Melted.level.order, y = value))+ geom_jitter(aes(color = Model), 
                                                                        position = position_jitter(0), size = 2.5) +
  xlab("") + ylab("Transition Rate")+ggtitle("Strict Standard Dataset: Pagel")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_colour_viridis_d(option = "viridis")+theme_classic()+coord_flip(ylim = c(0,0.03))

dev.size("in") #[1] 10.77778  5.87500

#Relaxed Standard

relaxed.Limbless.fit.fixed.dep.y.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Limbless.fit.fixed.dep.y.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Limbless.fit.fixed.dep.x.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Limbless.fit.fixed.dep.x.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Limbless.fit.fixed.dep.xy.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Limbless.fit.fixed.dep.xy.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Limbless.fit.fixed.dep.y.ER.Q <- as.data.frame(t(as.matrix(relaxed.Limbless.fit.fixed.dep.y.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Limbless.fit.fixed.dep.x.ER.Q <- as.data.frame(t(as.matrix(relaxed.Limbless.fit.fixed.dep.x.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Limbless.fit.fixed.dep.xy.ER.Q <- as.data.frame(t(as.matrix(relaxed.Limbless.fit.fixed.dep.xy.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))

#Export table
#relaxed.Limbless.fit.fixed.dep.y.ARD.Q <- setnames(relaxed.Limbless.fit.fixed.dep.y.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Limbless.fit.fixed.dep.x.ARD.Q <- setnames(relaxed.Limbless.fit.fixed.dep.x.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Limbless.fit.fixed.dep.xy.ARD.Q <- setnames(relaxed.Limbless.fit.fixed.dep.xy.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Limbless.fit.fixed.dep.y.ER.Q <- setnames(relaxed.Limbless.fit.fixed.dep.y.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Limbless.fit.fixed.dep.x.ER.Q <- setnames(relaxed.Limbless.fit.fixed.dep.x.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Limbless.fit.fixed.dep.xy.ER.Q <- setnames(relaxed.Limbless.fit.fixed.dep.xy.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#
#
#relaxed.Limbless.Q.export <- rbind(relaxed.Limbless.fit.fixed.dep.y.ARD.Q,
#                            relaxed.Limbless.fit.fixed.dep.x.ARD.Q,
#                            relaxed.Limbless.fit.fixed.dep.xy.ARD.Q,
#                            relaxed.Limbless.fit.fixed.dep.y.ER.Q,
#                            relaxed.Limbless.fit.fixed.dep.x.ER.Q,
#                            relaxed.Limbless.fit.fixed.dep.xy.ER.Q)
#
#relaxed.Limbless.Q.export <-as.data.frame(relaxed.Limbless.Q.export, row.names = c("ARD dependent y",
#                                                                     "ARD dependent x",
#                                                                     "ARD dependent x&y",
#                                                                     "ER dependent y",
#                                                                     "ER dependent x",
#                                                                     "ER dependent x&y"))
#
#relaxed.Limbless.AIC <- c(relaxed.Limbless.fit.fixed.dep.y.ARD$dependent.AIC,
#                          relaxed.Limbless.fit.fixed.dep.x.ARD$dependent.AIC,
#                          relaxed.Limbless.fit.fixed.dep.xy.ARD$dependent.AIC,
#                          relaxed.Limbless.fit.fixed.dep.y.ER$dependent.AIC,
#                          relaxed.Limbless.fit.fixed.dep.x.ER$dependent.AIC,
#                          relaxed.Limbless.fit.fixed.dep.xy.ER$dependent.AIC)
#
#
#
### Add in AIC scores
#relaxed.Limbless.Q.export$AIC<-relaxed.Limbless.AIC 
#
#relaxed.Limbless.Q.export

#####

relaxed.Limbless.fit.fixed.dep.y.ARD.Q <- setnames(relaxed.Limbless.fit.fixed.dep.y.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Limbless.fit.fixed.dep.x.ARD.Q <- setnames(relaxed.Limbless.fit.fixed.dep.x.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Limbless.fit.fixed.dep.xy.ARD.Q <- setnames(relaxed.Limbless.fit.fixed.dep.xy.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Limbless.fit.fixed.dep.y.ER.Q <- setnames(relaxed.Limbless.fit.fixed.dep.y.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Limbless.fit.fixed.dep.x.ER.Q <- setnames(relaxed.Limbless.fit.fixed.dep.x.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Limbless.fit.fixed.dep.xy.ER.Q <- setnames(relaxed.Limbless.fit.fixed.dep.xy.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))

relaxed.Limbless.Q <- rbind(relaxed.Limbless.fit.fixed.dep.y.ARD.Q,
                           relaxed.Limbless.fit.fixed.dep.x.ARD.Q,
                           relaxed.Limbless.fit.fixed.dep.xy.ARD.Q,
                           relaxed.Limbless.fit.fixed.dep.y.ER.Q,
                           relaxed.Limbless.fit.fixed.dep.x.ER.Q,
                           relaxed.Limbless.fit.fixed.dep.xy.ER.Q)

relaxed.Limbless.Q <-as.data.frame(relaxed.Limbless.Q, row.names = c("ARD dependent y",
                                                                   "ARD dependent x",
                                                                   "ARD dependent x&y",
                                                                   "ER dependent y",
                                                                   "ER dependent x",
                                                                   "ER dependent x&y"))

## Finalize table
setDT(relaxed.Limbless.Q, keep.rownames = TRUE)
colnames(relaxed.Limbless.Q)[1] <- "Model"

library(reshape2)
relaxed.Limbless.Q.Melted <- melt((relaxed.Limbless.Q))

relaxed.Limbless.Q.Melted.level.order <- factor(relaxed.Limbless.Q.Melted$variable, level = c('PadlessNonArboreal.PadbearingNonArboreal',
                                                                                              'PadlessArboreal.PadbearingArboreal',
                                                                                              'PadbearingArboreal.PadlessArboreal',
                                                                                              'PadbearingNonArboreal.PadlessNonArboreal',
                                                                                              'PadbearingArboreal.PadbearingNonArboreal',
                                                                                              'PadlessNonArboreal.PadlessArboreal',
                                                                                              'PadlessArboreal.PadlessNonArboreal', 
                                                                                              'PadbearingNonArboreal.PadbearingArboreal'))


Relaxed.Standard.Q.Jitter <- ggplot(relaxed.Limbless.Q.Melted, aes(x = relaxed.Limbless.Q.Melted.level.order, y = value))+ geom_jitter(aes(color = Model), 
                                                                            position = position_jitter(0), size = 2.5) +
  xlab("") + ylab("Transition Rate")+ggtitle("Relaxed Standard Dataset: Pagel")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_colour_viridis_d(option = "viridis")+theme_classic()+coord_flip(ylim = c(0,0.03))

dev.size("in") #[1] 10.77778  5.87500


## Together all
ggarrange(
  Strict.Standard.Q.Jitter,Relaxed.Standard.Q.Jitter, Strict.Full.Q.Jitter, Relaxed.Full.Q.Jitter, labels = c("A", "B", "C", "D"),
  common.legend = TRUE, legend = "top")
together + rremove("ggtitle")

## Only best-fitting models from the four dataset

Best.Models.Pagel.Q <- rbind(strict.Limbless.fit.fixed.dep.xy.ARD.Q,
                            relaxed.Limbless.fit.fixed.dep.y.ARD.Q,
                            strict.Full.fit.fixed.dep.xy.ARD.Q,
                            relaxed.Full.fit.fixed.dep.xy.ARD.Q)

Best.Models.Pagel.Q <-as.data.frame(Best.Models.Pagel.Q, row.names = c("Standard Strict",
                                                                     "Standard Relaxed",
                                                                     "Full Strict",
                                                                     "Full Relaxed"))
#Add p-values and AIC scores 
Best.Models.Pagel.Q.export <- Best.Models.Pagel.Q
Best.Models.Pagel.Q.export$p.value <-  rbind(strict.Limbless.fit.fixed.dep.xy.ARD$P[1],
                                      relaxed.Limbless.fit.fixed.dep.y.ARD$P[1],
                                      strict.Full.fit.fixed.dep.xy.ARD$P[1],
                                      relaxed.Full.fit.fixed.dep.xy.ARD$P[1])

Best.Models.Pagel.Q.export$log.lik <-  rbind(strict.Limbless.fit.fixed.dep.xy.ARD$dependent.logL[1],
                                             relaxed.Limbless.fit.fixed.dep.y.ARD$dependent.logL[1],
                                             strict.Full.fit.fixed.dep.xy.ARD$dependent.logL[1],
                                             relaxed.Full.fit.fixed.dep.xy.ARD$dependent.logL[1])
write.csv(Best.Models.Pagel.Q.export, "Best.Models.Pagel.SuppMat.Rates.P.14Feb21.csv")

## Finalize table
setDT(Best.Models.Pagel.Q, keep.rownames = TRUE)
colnames(Best.Models.Pagel.Q)[1] <- "Dataset"
library(reshape2)
Best.Models.Pagel.Q.Melted <- melt((Best.Models.Pagel.Q))

Best.Models.Pagel.Q.Melted.Jitter <- ggplot(Best.Models.Pagel.Q.Melted, aes(x = reorder(variable,value), y = value))+ geom_jitter(aes(color = Dataset), 
                                                                                                          position = position_jitter(0.05), size = 3) +
  xlab("") + ylab("Transition Rate")+ggtitle("Pagel AICw Winners")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_color_manual(values=c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"))+theme_classic()+coord_flip()
dev.size("in") #[1] 10.77778  5.87500

Best.Models.MuHiSSE.Q.Melted.Jitter <- ggplot(Best.Models.MuHiSSE.Q.Melted, aes(x=SSE.Best.Models.level.order, y = value))+ geom_jitter(aes(color = Dataset), 
                                                                                                                                        position = position_jitter(0.05), size = 3) +
  xlab("") + ylab("Transition Rate")+ggtitle("MuHiSSE AICw Winners")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_color_manual(values=c( "#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"))+theme_classic()+coord_flip()


## Not included-- multiplying rates by the sum of edge lengths to obtain a rough estimate of the number of state transitions

#Pagel converted by branch lengths
#Best.Models.Pagel.Q.Parameter.Estimates <- Best.Models.Pagel.Q.Parameter %>% column_to_rownames(var="Dataset")
#Best.Models.Pagel.Q.Parameter.Estimates.Standard <- Best.Models.Pagel.Q.Parameter.Estimates[1:2,]*sum(LimblessTree$edge.length)
#Best.Models.Pagel.Q.Parameter.Estimates.Full <- Best.Models.Pagel.Q.Parameter.Estimates[3:4,]*sum(FullTree$edge.length)
#Parameters.Best.Models.Pagel.Q.Combined <- rbind(Best.Models.Pagel.Q.Parameter.Estimates.Standard,
#                                                 Best.Models.Pagel.Q.Parameter.Estimates.Full)
#Parameters.Best.Models.Pagel.Q.Combined <- tibble::rownames_to_column(Parameters.Best.Models.Pagel.Q.Combined, "Dataset")
#Parameters.Best.Models.Pagel.Q.Combined.Melted <- melt((Parameters.Best.Models.Pagel.Q.Combined))
#
#Parameters.Best.Models.Pagel.Q.Combined.Melted.Jitter <- ggplot(Parameters.Best.Models.Pagel.Q.Combined.Melted, aes(x = reorder(variable,value), y = value))+ geom_jitter(aes(color = Dataset), 
#                                                                                                                                  position = position_jitter(0.2), size = 3) +
#  xlab("") + ylab("Number of Transitions")+ggtitle("Pagel AICw Winners: Converted")+
#  stat_summary(aes(color = value), size = 0.3,
#               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_color_manual(values=c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"))+theme_classic()+coord_flip()
#
##MuHiSSE converted by branch lengths
#Best.Models.MuHiSSE <- tibble::rownames_to_column(Best.Models.MuHiSSE, "Dataset")
#Best.Models.MuHiSSE.Q.Parameter <- Best.Models.MuHiSSE
#Best.Models.MuHiSSE.Q.Parameter.Estimates <- Best.Models.MuHiSSE.Q.Parameter %>% column_to_rownames(var="Dataset")
#Best.Models.MuHiSSE.Q.Parameter.Estimates.Standard <- Best.Models.MuHiSSE.Q.Parameter.Estimates[1:2,]*sum(LimblessTree$edge.length)
#Best.Models.MuHiSSE.Q.Parameter.Estimates.Full <- Best.Models.MuHiSSE.Q.Parameter.Estimates[3:4,]*sum(FullTree$edge.length)
#Parameters.Best.Models.MuHiSSE.Q.Combined <- rbind(Best.Models.MuHiSSE.Q.Parameter.Estimates.Standard,
#                                                 Best.Models.MuHiSSE.Q.Parameter.Estimates.Full)
#Parameters.Best.Models.MuHiSSE.Q.Combined <- tibble::rownames_to_column(Parameters.Best.Models.MuHiSSE.Q.Combined, "Dataset")
#Parameters.Best.Models.MuHiSSE.Q.Combined.Melted <- melt((Parameters.Best.Models.MuHiSSE.Q.Combined))
#
#Parameters.Best.Models.MuHiSSE.Q.Combined.Melted.level.order <- factor(Parameters.Best.Models.MuHiSSE.Q.Combined.Melted$variable, level = c('q00A_01A',
#                                                                                       'q10A_11A',
#                                                                                       'q11A_10A',
#                                                                                       'q01A_00A',
#                                                                                       'q00A_10A',
#                                                                                       'q11A_01A',
#                                                                                       'q01A_11A',
#                                                                                       'q10A_00A'))
#
#
#Parameters.Best.Models.MuHiSSE.Q.Combined.Melted.Jitter <- ggplot(Parameters.Best.Models.MuHiSSE.Q.Combined.Melted, aes(x=Parameters.Best.Models.MuHiSSE.Q.Combined.Melted.level.order, y =value))+ geom_jitter(aes(color = Dataset), 
#                                                                                                                                                                          position = position_jitter(0.2), size = 3) +
#  xlab("") + ylab("Number of Transitions")+ggtitle("MuHiSSE AICw Winners: Converted")+
#  stat_summary(aes(color = value), size = 0.3,
#               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_color_manual(values=c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"))+theme_classic()+coord_flip()
#
###Together converted parameters
ggarrange(Best.Models.MuHiSSE.Q.Melted.Jitter+coord_flip(ylim = c(0,0.025)),
          Best.Models.Pagel.Q.Melted.Jitter+coord_flip(ylim = c(0,0.025)),
         align = "v",labels = c("A", "B"), ncol = 1, nrow=2, common.legend = TRUE,legend = "top")
#
#ggarrange(Parameters.Best.Models.MuHiSSE.Q.Combined.Melted.Jitter+coord_flip(ylim = c(0,1000)),
#          Parameters.Best.Models.Pagel.Q.Combined.Melted.Jitter+coord_flip(ylim = c(0,1000)),
#          align = "v",labels = c("A", "B"), ncol = 1, nrow=2,legend = "right")


## Compare SSE to fitPagel rate parameters
compare <- read.csv("MuHiSSEvsPagel.csv")
compare_transpose <- as.data.frame(t(as.matrix(compare)))
var_q00._.10_All <- var(as.integer(compare_transpose[2,]))
var(as.integer(compare_transpose[2,]))
var(as.integer(compare_transpose[3,]))
var(as.integer(compare_transpose[4,]))
var(as.integer(compare_transpose[5,]))
var(as.integer(compare_transpose[6,]))

View(compare_transpose[2,])

aggregate(compare, list(compare$q00._.10), FUN = var)



### Aggregate all rates
relaxed.Limbless.Q.export
strict.Limbless.Q.export
strict.Full.Q.export
Relaxed.Full.Q.Export

write.csv(relaxed.Limbless.Q.export, "Standard_Relaxed_Rates_AllModels.csv")
write.csv(strict.Limbless.Q.export, "Standard_Strict_Rates_AllModels.csv")
write.csv(Relaxed.Full.Q.Export, "Full_Relaxed_Rates_AllModels.csv")
write.csv(strict.Full.Q.export, "Full_Strict_Rates_AllModels.csv")

AllRates_Pagel <- cbind(relaxed.Limbless.Q.export,
                    strict.Limbless.Q.export,
                    strict.Full.Q.export,
                    Relaxed.Full.Q.Export)




