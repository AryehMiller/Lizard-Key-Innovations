## fitPagel with fixed ("fossilized") root state
## Full and Standard datasets

library(phytools)
library(ape)
library(data.table)
library(ggplot2)
library(viridis)
library(Hmisc)
library(ggpubr)
library(geiger)
library(reshape2)

##################################
######                      ######
######   Strict "Standard"  ######
######                      ######
##################################

## Read in trait dataset

StandardDataset<- read.csv("StandardDataset_MeiriData.csv", row.names = 1)

## Read in phylogenies
StandardTree <- read.nexus(file = "StandardTree.tre")

  
## Add tip of zero-length to the root using "bind.tip" and then assigning it non-arboreal and padless
Standard.tree.root <-bind.tip(StandardTree,tip.label="ROOT",edge.length=0,
                          where=Ntip(StandardTree)+1)
#plotTree(Standard.tree.root,ftype="off",lwd=1)
nnL<-which(Standard.tree.root$tip.label=="ROOT")
#tiplabels(Standard.tree.root$tip.label[nnL],nnL,adj=c(-0.1,0.5),
#          frame="none",cex=0.8,font=3)

## Add fossilized root states
arboreal.strict.Standard<-as.factor(setNames(StandardDataset[,4],rownames(StandardDataset)))
toepads.Standard<-as.factor(setNames(StandardDataset[,7],rownames(StandardDataset)))

Standard.arboreal.strict.fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                  setNames(as.character(arboreal.strict.Standard),names(arboreal.strict.Standard))))
Standard.padless.fossilized <-as.factor(c(setNames("0","ROOT"),
                                 setNames(as.character(toepads.Standard),names(toepads.Standard))))

#fit fixed-root models for Strict Standard dataset

##ARD

## Dependent y model
strict.Standard.fit.fixed.dep.y.ARD<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.strict.fossilized,dep.var="y",model="ARD") #strict

## Dependent x model
strict.Standard.fit.fixed.dep.x.ARD <-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.strict.fossilized,dep.var="x",model="ARD") #strict

## Interdependent model
strict.Standard.fit.fixed.dep.xy.ARD<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.strict.fossilized,dep.var="xy",model="ARD") #strict


## ER

## Dependent y model
strict.Standard.fit.fixed.dep.y.ER<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.strict.fossilized,dep.var="y",model="ER") #strict

## Dependent x model
strict.Standard.fit.fixed.dep.x.ER <-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.strict.fossilized,dep.var="x",model="ER") #strict

## Interdependent model
strict.Standard.fit.fixed.dep.xy.ER<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.strict.fossilized,dep.var="xy",model="ER") #strict

strict.Standard.aic <-setNames(c(strict.Standard.fit.fixed.dep.y.ARD$dependent.AIC,
                                  strict.Standard.fit.fixed.dep.x.ARD$dependent.AIC,
                                  strict.Standard.fit.fixed.dep.xy.ARD$dependent.AIC,
                                  strict.Standard.fit.fixed.dep.y.ER$dependent.AIC,
                                  strict.Standard.fit.fixed.dep.x.ER$dependent.AIC,
                                  strict.Standard.fit.fixed.dep.xy.ER$dependent.AIC),
                                c("ARD dependent y",
                                  "ARD dependent x",
                                  "ARD dependent x&y",
                                  "ER dependent y",
                                  "ER dependent x",
                                  "ER dependent x&y"))

print(strict.Standard.aic)
strict.Standard.aicw <- aicw(strict.Standard.aic)
as.data.frame(strict.Standard.aicw$w)
plot(strict.Standard.fit.fixed.dep.xy.ARD)
##################################
######                      ######
######  Relaxed "Standard"  ######
######                      ######
##################################

## Add fossilized root states
arboreal.relaxed.Standard<-as.factor(setNames(StandardDataset[,5],rownames(StandardDataset)))

Standard.arboreal.relaxed.fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                                  setNames(as.character(arboreal.relaxed.Standard),names(arboreal.relaxed.Standard))))

## Fit fixed-root models for Standard Relaxed dataset

## ARD

## Dependent y model
relaxed.Standard.fit.fixed.dep.y.ARD<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.relaxed.fossilized,dep.var="y",model="ARD") #strict

#relaxed.Standard.fit.fixed.dep.y.ARD
#plot(relaxed.Standard.fit.fixed.dep.y.ARD)
#title("relaxed.Standard.fit.fixed.dep.y.ARD")

## Dependent x model
relaxed.Standard.fit.fixed.dep.x.ARD <-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.relaxed.fossilized,dep.var="x",model="ARD") #strict

## Interdependent model
relaxed.Standard.fit.fixed.dep.xy.ARD<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.relaxed.fossilized,dep.var="xy",model="ARD") #strict

## ER

## Dependent y model
relaxed.Standard.fit.fixed.dep.y.ER<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.relaxed.fossilized,dep.var="y",model="ER") #strict

## Dependent x model
relaxed.Standard.fit.fixed.dep.x.ER <-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.relaxed.fossilized,dep.var="x",model="ER") #strict

## Interdependent model
relaxed.Standard.fit.fixed.dep.xy.ER<-fitPagel(Standard.tree.root,Standard.padless.fossilized,Standard.arboreal.relaxed.fossilized,dep.var="xy",model="ER") #strict

## Organize AIC scores

relaxed.Standard.aic <-setNames(c(relaxed.Standard.fit.fixed.dep.y.ARD$dependent.AIC,
                                     relaxed.Standard.fit.fixed.dep.x.ARD$dependent.AIC,
                                     relaxed.Standard.fit.fixed.dep.xy.ARD$dependent.AIC,
                                     relaxed.Standard.fit.fixed.dep.y.ER$dependent.AIC,
                                     relaxed.Standard.fit.fixed.dep.x.ER$dependent.AIC,
                                     relaxed.Standard.fit.fixed.dep.xy.ER$dependent.AIC),
                                   c("ARD dependent y",
                                     "ARD dependent x",
                                     "ARD dependent x&y",
                                     "ER dependent y",
                                     "ER dependent x",
                                     "ER dependent x&y"))

print(relaxed.Standard.aic)
aicw(relaxed.Standard.aic)
relaxed.Standard.aicw <- aicw(relaxed.Standard.aic)
as.data.frame(relaxed.Standard.aicw$w)
plot(relaxed.Standard.fit.fixed.dep.y.ARD)
##################################
######                      ######
######   Strict "Full"      ######
######                      ######
##################################

#Full trait dataset
FullDataset<- read.csv("FullDataset_MeiriData.csv", row.names = 1)

#Full phylogeny
FullTree <- read.nexus(file = "FullTree.tre")

## Add tip of zero-length to the root using "bind.tip" and then assigning it non-arboreal and padless
Full.tree.root <-bind.tip(FullTree,tip.label="ROOT",edge.length=0,
                              where=Ntip(FullTree)+1)
#plotTree(Full.tree.root,ftype="off",lwd=1)
nnF<-which(Full.tree.root$tip.label=="ROOT")
#tiplabels(Full.tree.root$tip.label[nnF],nnF,adj=c(-0.1,0.5),
#          frame="none",cex=0.8,font=3)

## Add fossilized root states
arboreal.strict.Full<-as.factor(setNames(FullDataset[,4],rownames(FullDataset)))
toepads.Full<-as.factor(setNames(FullDataset[,7],rownames(FullDataset)))

Full.arboreal.strict.fossilized <-as.factor(c(setNames("not.arboreal","ROOT"),
                                                  setNames(as.character(arboreal.strict.Full),names(arboreal.strict.Full))))
Full.padless.fossilized <-as.factor(c(setNames("0","ROOT"),
                                          setNames(as.character(toepads.Full),names(toepads.Full))))

#fit fixed-root models for Strict Standard dataset

##ARD

## Dependent y model
strict.Full.fit.fixed.dep.y.ARD<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="y",model="ARD") #strict

## Dependent x model
strict.Full.fit.fixed.dep.x.ARD <-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="x",model="ARD") #strict

## Interdependent model
strict.Full.fit.fixed.dep.xy.ARD<-fitPagel(Full.tree.root,Full.padless.fossilized,Full.arboreal.strict.fossilized,dep.var="xy",model="ARD") #strict

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
strict.Full.aicw <- aicw(strict.Full.aic)
as.data.frame(strict.Full.aicw$w)
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
relaxed.Full.aicw <- aicw(relaxed.Full.aic)
as.data.frame(relaxed.Full.aicw$w)
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

strict.Standard.fit.fixed.dep.y.ARD.Q <- as.data.frame(t(as.matrix(strict.Standard.fit.fixed.dep.y.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Standard.fit.fixed.dep.x.ARD.Q <- as.data.frame(t(as.matrix(strict.Standard.fit.fixed.dep.x.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Standard.fit.fixed.dep.xy.ARD.Q <- as.data.frame(t(as.matrix(strict.Standard.fit.fixed.dep.xy.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Standard.fit.fixed.dep.y.ER.Q <- as.data.frame(t(as.matrix(strict.Standard.fit.fixed.dep.y.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Standard.fit.fixed.dep.x.ER.Q <- as.data.frame(t(as.matrix(strict.Standard.fit.fixed.dep.x.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
strict.Standard.fit.fixed.dep.xy.ER.Q <- as.data.frame(t(as.matrix(strict.Standard.fit.fixed.dep.xy.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))

### Table for export

#strict.Standard.fit.fixed.dep.y.ARD.Q <- setnames(strict.Standard.fit.fixed.dep.y.ARD.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Standard.fit.fixed.dep.x.ARD.Q <- setnames(strict.Standard.fit.fixed.dep.x.ARD.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Standard.fit.fixed.dep.xy.ARD.Q <- setnames(strict.Standard.fit.fixed.dep.xy.ARD.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Standard.fit.fixed.dep.y.ER.Q <- setnames(strict.Standard.fit.fixed.dep.y.ER.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Standard.fit.fixed.dep.x.ER.Q <- setnames(strict.Standard.fit.fixed.dep.x.ER.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#strict.Standard.fit.fixed.dep.xy.ER.Q <- setnames(strict.Standard.fit.fixed.dep.xy.ER.Q,c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#
#strict.Standard.Q.export <- rbind(strict.Standard.fit.fixed.dep.y.ARD.Q,
#                           strict.Standard.fit.fixed.dep.x.ARD.Q,
#                           strict.Standard.fit.fixed.dep.xy.ARD.Q,
#                           strict.Standard.fit.fixed.dep.y.ER.Q,
#                           strict.Standard.fit.fixed.dep.x.ER.Q,
#                           strict.Standard.fit.fixed.dep.xy.ER.Q)
#
#strict.Standard.Q.export <-as.data.frame(strict.Standard.Q.export, row.names = c("ARD dependent y",
#                                                                   "ARD dependent x",
#                                                                   "ARD dependent x&y",
#                                                                   "ER dependent y",
#                                                                   "ER dependent x",
#                                                                   "ER dependent x&y"))
#
#strict.Standard.AIC <- c(strict.Standard.fit.fixed.dep.y.ARD$dependent.AIC,
#                         strict.Standard.fit.fixed.dep.x.ARD$dependent.AIC,
#                         strict.Standard.fit.fixed.dep.xy.ARD$dependent.AIC,
#                         strict.Standard.fit.fixed.dep.y.ER$dependent.AIC,
#                         strict.Standard.fit.fixed.dep.x.ER$dependent.AIC,
#                         strict.Standard.fit.fixed.dep.xy.ER$dependent.AIC)
#
### Add in AIC scores
#strict.Standard.Q.export$AIC<-strict.Standard.AIC 
#
#strict.Standard.Q.export

###

strict.Standard.fit.fixed.dep.y.ARD.Q <- setnames(strict.Standard.fit.fixed.dep.y.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Standard.fit.fixed.dep.x.ARD.Q <- setnames(strict.Standard.fit.fixed.dep.x.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Standard.fit.fixed.dep.xy.ARD.Q <- setnames(strict.Standard.fit.fixed.dep.xy.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Standard.fit.fixed.dep.y.ER.Q <- setnames(strict.Standard.fit.fixed.dep.y.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Standard.fit.fixed.dep.x.ER.Q <- setnames(strict.Standard.fit.fixed.dep.x.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
strict.Standard.fit.fixed.dep.xy.ER.Q <- setnames(strict.Standard.fit.fixed.dep.xy.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))

strict.Standard.Q <- rbind(strict.Standard.fit.fixed.dep.y.ARD.Q,
                       strict.Standard.fit.fixed.dep.x.ARD.Q,
                       strict.Standard.fit.fixed.dep.xy.ARD.Q,
                       strict.Standard.fit.fixed.dep.y.ER.Q,
                       strict.Standard.fit.fixed.dep.x.ER.Q,
                       strict.Standard.fit.fixed.dep.xy.ER.Q)

strict.Standard.Q <-as.data.frame(strict.Standard.Q, row.names = c("ARD dependent y",
                                                           "ARD dependent x",
                                                           "ARD dependent x&y",
                                                           "ER dependent y",
                                                           "ER dependent x",
                                                           "ER dependent x&y"))



## Finalize table
setDT(strict.Standard.Q, keep.rownames = TRUE)
colnames(strict.Standard.Q)[1] <- "Model"

strict.Standard.Q.Melted <- melt((strict.Standard.Q))

strict.Standard.Q.Melted.level.order <- factor(strict.Standard.Q.Melted$variable, level = c('PadlessNonArboreal.PadbearingNonArboreal',
                                                                                              'PadlessArboreal.PadbearingArboreal',
                                                                                              'PadbearingArboreal.PadlessArboreal',
                                                                                              'PadbearingNonArboreal.PadlessNonArboreal',
                                                                                              'PadbearingArboreal.PadbearingNonArboreal',
                                                                                              'PadlessNonArboreal.PadlessArboreal',
                                                                                              'PadlessArboreal.PadlessNonArboreal', 
                                                                                              'PadbearingNonArboreal.PadbearingArboreal'))

Strict.Standard.Q.Jitter <-ggplot(strict.Standard.Q.Melted, aes(x = strict.Standard.Q.Melted.level.order, y = value))+ geom_jitter(aes(color = Model), 
                                                                        position = position_jitter(0), size = 2.5) +
  xlab("") + ylab("Transition Rate")+ggtitle("Strict Standard Dataset: Pagel")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_colour_viridis_d(option = "viridis")+theme_classic()+coord_flip(ylim = c(0,0.03))

dev.size("in") #[1] 10.77778  5.87500

#Relaxed Standard

relaxed.Standard.fit.fixed.dep.y.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Standard.fit.fixed.dep.y.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Standard.fit.fixed.dep.x.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Standard.fit.fixed.dep.x.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Standard.fit.fixed.dep.xy.ARD.Q <- as.data.frame(t(as.matrix(relaxed.Standard.fit.fixed.dep.xy.ARD$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Standard.fit.fixed.dep.y.ER.Q <- as.data.frame(t(as.matrix(relaxed.Standard.fit.fixed.dep.y.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Standard.fit.fixed.dep.x.ER.Q <- as.data.frame(t(as.matrix(relaxed.Standard.fit.fixed.dep.x.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))
relaxed.Standard.fit.fixed.dep.xy.ER.Q <- as.data.frame(t(as.matrix(relaxed.Standard.fit.fixed.dep.xy.ER$dependent.Q[c(2,3,5,8,9,12,14,15)])))

#Export table
#relaxed.Standard.fit.fixed.dep.y.ARD.Q <- setnames(relaxed.Standard.fit.fixed.dep.y.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Standard.fit.fixed.dep.x.ARD.Q <- setnames(relaxed.Standard.fit.fixed.dep.x.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Standard.fit.fixed.dep.xy.ARD.Q <- setnames(relaxed.Standard.fit.fixed.dep.xy.ARD.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Standard.fit.fixed.dep.y.ER.Q <- setnames(relaxed.Standard.fit.fixed.dep.y.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Standard.fit.fixed.dep.x.ER.Q <- setnames(relaxed.Standard.fit.fixed.dep.x.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#relaxed.Standard.fit.fixed.dep.xy.ER.Q <- setnames(relaxed.Standard.fit.fixed.dep.xy.ER.Q, c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
#
#
#relaxed.Standard.Q.export <- rbind(relaxed.Standard.fit.fixed.dep.y.ARD.Q,
#                            relaxed.Standard.fit.fixed.dep.x.ARD.Q,
#                            relaxed.Standard.fit.fixed.dep.xy.ARD.Q,
#                            relaxed.Standard.fit.fixed.dep.y.ER.Q,
#                            relaxed.Standard.fit.fixed.dep.x.ER.Q,
#                            relaxed.Standard.fit.fixed.dep.xy.ER.Q)
#
#relaxed.Standard.Q.export <-as.data.frame(relaxed.Standard.Q.export, row.names = c("ARD dependent y",
#                                                                     "ARD dependent x",
#                                                                     "ARD dependent x&y",
#                                                                     "ER dependent y",
#                                                                     "ER dependent x",
#                                                                     "ER dependent x&y"))
#
#relaxed.Standard.AIC <- c(relaxed.Standard.fit.fixed.dep.y.ARD$dependent.AIC,
#                          relaxed.Standard.fit.fixed.dep.x.ARD$dependent.AIC,
#                          relaxed.Standard.fit.fixed.dep.xy.ARD$dependent.AIC,
#                          relaxed.Standard.fit.fixed.dep.y.ER$dependent.AIC,
#                          relaxed.Standard.fit.fixed.dep.x.ER$dependent.AIC,
#                          relaxed.Standard.fit.fixed.dep.xy.ER$dependent.AIC)
#
#
#
### Add in AIC scores
#relaxed.Standard.Q.export$AIC<-relaxed.Standard.AIC 
#
#relaxed.Standard.Q.export

#####

relaxed.Standard.fit.fixed.dep.y.ARD.Q <- setnames(relaxed.Standard.fit.fixed.dep.y.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Standard.fit.fixed.dep.x.ARD.Q <- setnames(relaxed.Standard.fit.fixed.dep.x.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Standard.fit.fixed.dep.xy.ARD.Q <- setnames(relaxed.Standard.fit.fixed.dep.xy.ARD.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Standard.fit.fixed.dep.y.ER.Q <- setnames(relaxed.Standard.fit.fixed.dep.y.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Standard.fit.fixed.dep.x.ER.Q <- setnames(relaxed.Standard.fit.fixed.dep.x.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))
relaxed.Standard.fit.fixed.dep.xy.ER.Q <- setnames(relaxed.Standard.fit.fixed.dep.xy.ER.Q, old = c('V1','V2','V3','V4','V5','V6','V7','V8'), new = c("PadlessNonArboreal.PadlessArboreal","PadbearingArboreal.PadlessArboreal", "PadlessArboreal.PadlessNonArboreal","PadbearingNonArboreal.PadlessNonArboreal","PadlessArboreal.PadbearingArboreal","PadbearingNonArboreal.PadbearingArboreal","PadlessNonArboreal.PadbearingNonArboreal","PadbearingArboreal.PadbearingNonArboreal"))

relaxed.Standard.Q <- rbind(relaxed.Standard.fit.fixed.dep.y.ARD.Q,
                           relaxed.Standard.fit.fixed.dep.x.ARD.Q,
                           relaxed.Standard.fit.fixed.dep.xy.ARD.Q,
                           relaxed.Standard.fit.fixed.dep.y.ER.Q,
                           relaxed.Standard.fit.fixed.dep.x.ER.Q,
                           relaxed.Standard.fit.fixed.dep.xy.ER.Q)

relaxed.Standard.Q <-as.data.frame(relaxed.Standard.Q, row.names = c("ARD dependent y",
                                                                   "ARD dependent x",
                                                                   "ARD dependent x&y",
                                                                   "ER dependent y",
                                                                   "ER dependent x",
                                                                   "ER dependent x&y"))

## Finalize table
setDT(relaxed.Standard.Q, keep.rownames = TRUE)
colnames(relaxed.Standard.Q)[1] <- "Model"

relaxed.Standard.Q.Melted <- melt((relaxed.Standard.Q))

relaxed.Standard.Q.Melted.level.order <- factor(relaxed.Standard.Q.Melted$variable, level = c('PadlessNonArboreal.PadbearingNonArboreal',
                                                                                              'PadlessArboreal.PadbearingArboreal',
                                                                                              'PadbearingArboreal.PadlessArboreal',
                                                                                              'PadbearingNonArboreal.PadlessNonArboreal',
                                                                                              'PadbearingArboreal.PadbearingNonArboreal',
                                                                                              'PadlessNonArboreal.PadlessArboreal',
                                                                                              'PadlessArboreal.PadlessNonArboreal', 
                                                                                              'PadbearingNonArboreal.PadbearingArboreal'))


Relaxed.Standard.Q.Jitter <- ggplot(relaxed.Standard.Q.Melted, aes(x = relaxed.Standard.Q.Melted.level.order, y = value))+ geom_jitter(aes(color = Model), 
                                                                            position = position_jitter(0), size = 2.5) +
  xlab("") + ylab("Transition Rate")+ggtitle("Relaxed Standard Dataset: Pagel")+
  stat_summary(aes(color = value), size = 0.3,
               fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))+ scale_colour_viridis_d(option = "viridis")+theme_classic()+coord_flip(ylim = c(0,0.03))

dev.size("in") #[1] 10.77778  5.87500


## All together
ggarrange(
  Strict.Standard.Q.Jitter,Relaxed.Standard.Q.Jitter, Strict.Full.Q.Jitter, Relaxed.Full.Q.Jitter, labels = c("A", "B", "C", "D"),
  common.legend = TRUE, legend = "top")
together + rremove("ggtitle")

## Only best-fitting models from the four dataset

Best.Models.Pagel.Q <- rbind(strict.Standard.fit.fixed.dep.xy.ARD.Q,
                            relaxed.Standard.fit.fixed.dep.y.ARD.Q,
                            strict.Full.fit.fixed.dep.xy.ARD.Q,
                            relaxed.Full.fit.fixed.dep.xy.ARD.Q)

Best.Models.Pagel.Q <-as.data.frame(Best.Models.Pagel.Q, row.names = c("Standard Strict",
                                                                     "Standard Relaxed",
                                                                     "Full Strict",
                                                                     "Full Relaxed"))
#Add p-values and AIC scores for table
Best.Models.Pagel.Q.export <- Best.Models.Pagel.Q
Best.Models.Pagel.Q.export$p.value <-  rbind(strict.Standard.fit.fixed.dep.xy.ARD$P[1],
                                      relaxed.Standard.fit.fixed.dep.y.ARD$P[1],
                                      strict.Full.fit.fixed.dep.xy.ARD$P[1],
                                      relaxed.Full.fit.fixed.dep.xy.ARD$P[1])

Best.Models.Pagel.Q.export$log.lik <-  rbind(strict.Standard.fit.fixed.dep.xy.ARD$dependent.logL[1],
                                             relaxed.Standard.fit.fixed.dep.y.ARD$dependent.logL[1],
                                             strict.Full.fit.fixed.dep.xy.ARD$dependent.logL[1],
                                             relaxed.Full.fit.fixed.dep.xy.ARD$dependent.logL[1])

#write.csv(Best.Models.Pagel.Q.export, "Best.Models.Pagel.SuppMat.Rates.P.14Feb21.csv")

## Finalize table
setDT(Best.Models.Pagel.Q, keep.rownames = TRUE)
colnames(Best.Models.Pagel.Q)[1] <- "Dataset"

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

#### All datasets/models/rates

write.csv(relaxed.Standard.Q.export, "Standard_Relaxed_Rates_AllModels.csv")
write.csv(strict.Standard.Q.export, "Standard_Strict_Rates_AllModels.csv")
write.csv(Relaxed.Full.Q.Export, "Full_Relaxed_Rates_AllModels.csv")
write.csv(strict.Full.Q.export, "Full_Strict_Rates_AllModels.csv")
