###Getting the top most abundant OTUs for each host species
##For Matsuda et al. 2021

#load libraries
library(phylsoeq)

setwd("~/Documents/OSUDocs/Projects/ACE21/Hawaiian-coral-thermal-tolerance/Data/16s_analysis/")
load("16s_phyloseq4HE.RData")

bac.exp <- subset_samples(Bac.seq, Type == "sample")
  
##1. Montipora capitata
mcap <- subset_samples(Bac.seq, Species == "Montipora_capitata")
mcap.ra <- transform_sample_counts(mcap, function(x) x/ sum(x))

# Selecting top 10 OTUs
TopNOTUs = names(sort(taxa_sums(mcap.ra), TRUE)[1:10])
mcap.abund = prune_taxa(TopNOTUs, mcap.ra)
mcap.abund.melt <- psmelt(mcap.abund)
mcap.abund.melt$top.group <- "abund"
write.csv(mcap.abund.melt, "top_abundance/output/mcap_abund.csv")


##2. Porites compressa
pcomp <- subset_samples(Bac.seq, Species == "Porites_compressa")
pcomp.ra <- transform_sample_counts(pcomp, function(x) x/ sum(x))

# Selecting top 10 OTUs
TopNOTUs = names(sort(taxa_sums(pcomp.ra), TRUE)[1:10])
pcomp.abund = prune_taxa(TopNOTUs, pcomp.ra)
pcomp.abund.melt <- psmelt(pcomp.abund)
pcomp.abund.melt$top.group <- "abund"
write.csv(pcomp.abund.melt, "top_abundance/output/pcomp_abund.csv")


##3. Pavona varians
pvar <- subset_samples(Bac.seq, Species == "Pavona_varians")
pvar.ra <- transform_sample_counts(pvar, function(x) x/ sum(x))

# Selecting top 10 OTUs
TopNOTUs = names(sort(taxa_sums(pvar.ra), TRUE)[1:10])
pvar.abund = prune_taxa(TopNOTUs, pvar.ra)
pvar.abund.melt <- psmelt(pvar.abund)
pvar.abund.melt$top.group <- "abund"
write.csv(pvar.abund.melt, "top_abundance/output/pvar_abund.csv")


##4. Pocillopora acuta
pacu <- subset_samples(Bac.seq, Species == "Pocillopora_acuta")
levels(sample_data(pacu)$Time.Point)[levels(sample_data(pacu)$Time.Point)=="F1"] <- "T1"
pacu.ra <- transform_sample_counts(pacu, function(x) x/ sum(x))

# Selecting top 10 OTUs
TopNOTUs = names(sort(taxa_sums(pacu.ra), TRUE)[1:10])
pacu.abund = prune_taxa(TopNOTUs, pacu.ra)
pacu.abund.melt <- psmelt(pacu.abund)
pacu.abund.melt$top.group <- "abund"
write.csv(pacu.abund.melt, "top_abundance/output/pacu_abund.csv")


