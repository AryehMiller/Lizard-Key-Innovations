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

sse_dataset_strict_toepads <- read.csv("sse_dataset_strict_toepads.csv")
sse_dataset_relaxed_toepads <- read.csv("sse_dataset_relaxed_toepads.csv")

## Read in phylogeny

FullTree <- read.nexus(file = "pruned.tree.Meiri.RevisedPyron.tre")

## Full D-Stat Toepads

Concatenated.Full <- sse_dataset_strict_toepads$Concatenated

## Standard D-Stat Toepads

Toepads.Full <- sse_dataset_strict_toepads$toepads
Full.Toepads.Combo <- data.frame(Concatenated.Full, Toepads.Full)
Full.Toepads.PhyloD <- phylo.d(Full.Toepads.Combo,FullTree, names.col = Concatenated.Full, binvar = Toepads.Full)
Full.Toepads.PhyloD

## Full Strict Arboreality D-Stat

Full.Arboreal.Strict<- sse_dataset_strict_toepads$strict.arboreal
Full.Strict.Arboreal.Combo <- data.frame(Concatenated.Full, Full.Arboreal.Strict)
Full.Strict.Arboreal.PhyloD <- phylo.d(Full.Strict.Arboreal.Combo, FullTree, names.col = Concatenated.Full, binvar=Full.Arboreal.Strict)
Full.Strict.Arboreal.PhyloD

## Full Relaxed Arboreality D-Stat

Full.Arboreal.Relaxed <- sse_dataset_relaxed_toepads$relaxed.arboreal
Full.Relaxed.Arboreal.Combo <- data.frame(Concatenated.Full, Full.Arboreal.Relaxed)
Full.Relaxed.Arboreal.PhyloD <- phylo.d(Full.Relaxed.Arboreal.Combo, FullTree, names.col = Concatenated.Full, binvar=Full.Arboreal.Relaxed)
Full.Relaxed.Arboreal.PhyloD

## Standard (limbless removed) datasets

## Read in dataset

sse_dataset_Standard_Strict_toepads<-read.csv("Limbless_sse_dataset_strict_toepads.csv")
sse_dataset_Standard_Relaxed_toepads<-read.csv("Limbless_sse_dataset_relaxed_toepads.csv")

## Read in tree

StandardTree <- read.nexus(file = "Meiri23Dec20_LimblessRemoved.tre")

## Isolate species names

Concatenated.Standard <- sse_dataset_Standard_Strict_toepads$Concatenated

## Standard D-Stat Toepads

Toepads.Standard <- sse_dataset_Standard_Strict_toepads$toepads
Toepads.Standard.Combo <- data.frame(Concatenated.Standard, Toepads.Standard)
Toepads.Standard.PhyloD <- phylo.d(Toepads.Standard.Combo,StandardTree, names.col = Concatenated.Standard, binvar = Toepads.Standard)
Toepads.Standard.PhyloD

## Standard Strict Arboreality D-Stat

Standard.Arboreal.Strict<- sse_dataset_Standard_Strict_toepads$strict.arboreal
Standard.Strict.Arboreal.Combo <- data.frame(Concatenated.Standard, Standard.Arboreal.Strict)
Standard.Strict.Arboreal.PhyloD <- phylo.d(Standard.Strict.Arboreal.Combo, StandardTree, names.col = Concatenated.Standard, binvar=Standard.Arboreal.Strict)
Standard.Strict.Arboreal.PhyloD

## Standard Relaxed Arboreality D-Stat

Standard.Arboreal.Relaxed <- sse_dataset_Standard_Relaxed_toepads$relaxed.arboreal
Standard.Relaxed.Arboreal.Combo <- data.frame(Concatenated.Standard, Standard.Arboreal.Relaxed)
Standard.Relaxed.Arboreal.PhyloD <- phylo.d(Standard.Relaxed.Arboreal.Combo, StandardTree, names.col = Concatenated.Standard, binvar=Standard.Arboreal.Relaxed)
Standard.Relaxed.Arboreal.PhyloD

