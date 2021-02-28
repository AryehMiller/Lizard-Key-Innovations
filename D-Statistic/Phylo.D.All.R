##########################################################################
######                                                              ######
######   	Tests for Phylogenetic Signal with the D-statistic        ######
######                                                              ######
##########################################################################

## Load packages 
library(ape)
library(caper)

## Full

## Read in datasets

Full_Strict_Dataset <- read.csv("FullStrict_TraitAnalysis.csv")
Full_Relaxed_Dataset <- read.csv("FullRelaxed_TraitAnalysis.csv")

## Read in phylogeny

FullTree <- read.nexus(file = "FullTree.tre")

## Full D-Stat Toepads

Concatenated.Full <- Full_Strict_Dataset$Concatenated

## Standard D-Stat Toepads

Toepads.Full <- Full_Strict_Dataset$toepads
Full.Toepads.Combo <- data.frame(Concatenated.Full, Toepads.Full)
Full.Toepads.PhyloD <- phylo.d(Full.Toepads.Combo,FullTree, names.col = Concatenated.Full, binvar = Toepads.Full)
Full.Toepads.PhyloD

## Full Strict Arboreality D-Stat

Full.Arboreal.Strict<- Full_Strict_Dataset$strict.arboreal
Full.Strict.Arboreal.Combo <- data.frame(Concatenated.Full, Full.Arboreal.Strict)
Full.Strict.Arboreal.PhyloD <- phylo.d(Full.Strict.Arboreal.Combo, FullTree, names.col = Concatenated.Full, binvar=Full.Arboreal.Strict)
Full.Strict.Arboreal.PhyloD

## Full Relaxed Arboreality D-Stat

Full.Arboreal.Relaxed <- Full_Relaxed_Dataset$relaxed.arboreal
Full.Relaxed.Arboreal.Combo <- data.frame(Concatenated.Full, Full.Arboreal.Relaxed)
Full.Relaxed.Arboreal.PhyloD <- phylo.d(Full.Relaxed.Arboreal.Combo, FullTree, names.col = Concatenated.Full, binvar=Full.Arboreal.Relaxed)
Full.Relaxed.Arboreal.PhyloD

## Standard (limbless removed) datasets

## Read in dataset

Standard_Strict_Dataset<-read.csv("StandardStrict_TraitAnalysis.csv")
Standard_Relaxed_Dataset<-read.csv("StandardRelaxed_TraitAnalysis.csv")

## Read in tree

StandardTree <- read.nexus(file = "StandardTree.tre")

## Isolate species names

Concatenated.Standard <- Standard_Strict_Dataset$Concatenated

## Standard D-Stat Toepads

Toepads.Standard <- Standard_Strict_Dataset$toepads
Toepads.Standard.Combo <- data.frame(Concatenated.Standard, Toepads.Standard)
Toepads.Standard.PhyloD <- phylo.d(Toepads.Standard.Combo,StandardTree, names.col = Concatenated.Standard, binvar = Toepads.Standard)
Toepads.Standard.PhyloD

## Standard Strict Arboreality D-Stat

Standard.Arboreal.Strict<- Standard_Strict_Dataset$strict.arboreal
Standard.Strict.Arboreal.Combo <- data.frame(Concatenated.Standard, Standard.Arboreal.Strict)
Standard.Strict.Arboreal.PhyloD <- phylo.d(Standard.Strict.Arboreal.Combo, StandardTree, names.col = Concatenated.Standard, binvar=Standard.Arboreal.Strict)
Standard.Strict.Arboreal.PhyloD

## Standard Relaxed Arboreality D-Stat

Standard.Arboreal.Relaxed <- Standard_Relaxed_Dataset$relaxed.arboreal
Standard.Relaxed.Arboreal.Combo <- data.frame(Concatenated.Standard, Standard.Arboreal.Relaxed)
Standard.Relaxed.Arboreal.PhyloD <- phylo.d(Standard.Relaxed.Arboreal.Combo, StandardTree, names.col = Concatenated.Standard, binvar=Standard.Arboreal.Relaxed)
Standard.Relaxed.Arboreal.PhyloD
