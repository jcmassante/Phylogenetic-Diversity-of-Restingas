# =======================================================================================================
# R code used in the following publication:
# =======================================================================================================
# Title:   ENVIRONMENT AND EVOLUTIONARY HISTORY DEPICT PHYLOGENETIC ALPHA AND BETA DIVERSITY
# IN THE ATLANTIC COASTAL WHITE-SAND WOODLANDS

# Authors: Jhonny C Massante & Pille Gerhold
#======================================================================================

# Analyses need information from the community data (comm), and the dataset (data.full).
# The full code was split into four parts:

# (1) Preprocessing, (2) Statistical analysis, (3) Figures, (4) plant families analysis .

# This is part 1, and it includes the following processes:

# Standardisation of species names
# Phylogenetic tree
# Indexes of phylogenetic alpha and beta diversity


#======================================================================================

# All analyses were carried out using R Project. Users will need to 
# create their own R Project or working directory.

#======================================================================================

# Load required packages
library(Taxonstand)#Names standardisation
library(V.PhyloMaker)#Phylogenetic tree
library(ape)#Phylogenetic analysis
library(picante)#Phylogenetic diversity
library(betapart)#Phylogenetic turnover
library(openxlsx)#Open .xlsx datasets
library(here)#Find files in the R project specific folder

#======================================================================================
# Standardisation of species names
# Input: data frame with original species names (from community data) according to the Taxonstand package requirements
# Output: data frame with standardised species names
# NOTE: it may be necessary to do manual standardisation of some specific taxa

to.standardise <- read.xlsx(here("Data", "to.standardise.xlsx"))# Input
names <- TPL(to.standardise$Full.name)# Output

# Save the output
write.xlsx(names, "splist.xlsx")


#======================================================================================

# Phylogenetic tree
# Input: data frame with species list (standardised names) according to V.Phylomaker package requirements
# Output: a time calibrated phylogenetic tree

# Read required species list
splist <- read.xlsx(here("Data", "splist.xlsx"))# Input

# Building tree using phylo.maker function from V.Phylomaker package
Phylogenetic.tree <- phylo.maker(sp.list = splist, tree = GBOTB.extended, nodes = nodes.info.1, scenarios = "S3")# Output

# Exclude tree fern species
# There was a gymnosperm species ("Podocarpus sellowii"), which has been excluded manually.
tree.fern <- c("Cyathea_microdonta" , "Cyathea_leucofolis", "Cyathea_delgadii", "Cyathea_corcovadensis", "Cyathea_atrovirens")
tree <- drop.tip(tree, tree.fern)


#======================================================================================

# Calculation of response variables (phylogenetic alpha and beta diversity)

# Calculate phylogenetic alpha diversity as the standardised effect size of the
# mean pairwise distance (MPD) and mean nearest taxon distance (MNTD) using the Picante package
# Input: matrix of community data (species x sites) and a phylogenetic tree


# Phylogenetic tree
tree <- read.tree(here("Data", "Phylogenetic.tree"))# Input

# Community data
comm <- read.xlsx(here("Data", "comm.xlsx"), rowNames = T)# Input
colnames(comm) <- gsub("\\.", "-", colnames(comm))# Replace symbols to prevent error

# Prune the phylogenetic tree
tree.pruned<- prune.sample(comm, tree)

# Sort the columns of community data to be in the same order as the tip labels of the phylogeny
comm<- comm[, tree.pruned$tip.label]

# Create species distance matrix from the pruned tree
dist<-cophenetic(tree)


### Phylogenetic alpha diversity ###

# Calculate standardised effect size of mean pairwise distance (SESmpd)
sesmpd<- ses.mpd(as.matrix(comm), (dist), null.model = "taxa.labels", runs = 999)# Output
#write.table(sesmpd, "sesmpd.csv", sep=",")# To include in the main dataset

# Calculate standardised effect size of mean nearest taxon distance (SESmntd)
sesmntd <- ses.mntd(as.matrix(comm), (dist), null.model = "taxa.labels", runs = 999)# Output
#write.table(sesmntd, "sesmntd.csv", sep=",")# To include in the main dataset

#======================================================================================

### Phylogenetic beta diversity ###

# Phylogenetic beta diversity metrics from Picante package.

# Comdist calculates the expected phylogenetic distance
# separating two individuals or taxa drawn randomly from different communities.
# Comdistnt calculates the average phylogenetic distance to the most similar 
# taxon or individual in the other community for taxa or individuals in two communities.

# Input: distance matrix calculated from the phylogenetic tree

# Deep phylogenetic levels
phylobeta.dpw <- comdist(comm, dist)
#save(phylobeta.dpw, file = "phylobeta.dpw.RData")# Save for further analysis

# Shallow phylogenetic levels
phylobeta.dnn <- comdistnt(comm, dist)
#save(phylobeta.dnn, file = "phylobeta.dnn.RData")# Save for further analysis


# Phylogenetic beta diversity from Betapart package
# phylo.beta.pair calculates dissimilarity matrix accounting for phylogenetic beta diversity,
# measured as Sorensen derived pairwise phylogenetic dissimilarity.

phylo.restinga <- phylo.beta.pair(comm, tree)# Total phylobeta (turnover + nestedness)
simpson.true.phylo.beta <- phylo.restinga$phylo.beta.sim# Only turnover (Simpson which is unaffected by species richness)
#save(simpson.true.phylo.beta, file = "simpson.true.phylo.beta.RData")# Save for further analysis
