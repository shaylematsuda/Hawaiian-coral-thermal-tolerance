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

##Select top 10 OTUs in T0 only
mcap.T0 <- subset_samples(mcap, Time.Point == "T0")
mcap.T0.ra <- transform_sample_counts(mcap.T0, function(x) x/ sum(x))
TopNOTUs = names(sort(taxa_sums(mcap.T0.ra), TRUE)[1:10])
mcap.T0.abund = prune_taxa(TopNOTUs, mcap.T0.ra)
mcap.T0.abund.melt <- psmelt(mcap.T0.abund)
mcap.T0.abund.melt$top.group <- "abund"
mcap.T0.abund.melt <- mcap.T0.abund.melt %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(mcap.T0.abund.melt, "../../top_abundance/output/mcap_T0_abund.csv")

##2. Porites compressa
pcomp <- subset_samples(Bac.seq, Species == "Porites_compressa")
pcomp.ra <- transform_sample_counts(pcomp, function(x) x/ sum(x))

# Selecting top 10 OTUs
TopNOTUs = names(sort(taxa_sums(pcomp.ra), TRUE)[1:10])
pcomp.abund = prune_taxa(TopNOTUs, pcomp.ra)
pcomp.abund.melt <- psmelt(pcomp.abund)
pcomp.abund.melt$top.group <- "abund"
write.csv(pcomp.abund.melt, "top_abundance/output/pcomp_abund.csv")


##T0 only
pcomp.T0 <- subset_samples(pcomp, Time.Point == "T0")
pcomp.T0.ra <- transform_sample_counts(pcomp.T0, function(x) x/ sum(x))
TopNOTUs = names(sort(taxa_sums(pcomp.T0.ra), TRUE)[1:10])
pcomp.T0.abund = prune_taxa(TopNOTUs, pcomp.T0.ra)
pcomp.T0.abund.melt <- psmelt(pcomp.T0.abund)
pcomp.T0.abund.melt$top.group <- "abund"
pcomp.T0.abund.melt <- pcomp.T0.abund.melt %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(pcomp.T0.abund.melt, "../../top_abundance/output/pcomp_T0_abund.csv")


##3. Pavona varians
pvar <- subset_samples(Bac.seq, Species == "Pavona_varians")
pvar.ra <- transform_sample_counts(pvar, function(x) x/ sum(x))

# Selecting top 10 OTUs
TopNOTUs = names(sort(taxa_sums(pvar.ra), TRUE)[1:10])
pvar.abund = prune_taxa(TopNOTUs, pvar.ra)
pvar.abund.melt <- psmelt(pvar.abund)
pvar.abund.melt$top.group <- "abund"
write.csv(pvar.abund.melt, "top_abundance/output/pvar_abund.csv")

#T0 only
pvar.T0 <- subset_samples(pvar, Time.Point == "T0")
pvar.T0.ra <- transform_sample_counts(pvar.T0, function(x) x/ sum(x))
TopNOTUs = names(sort(taxa_sums(pvar.T0.ra), TRUE)[1:10])
pvar.T0.abund = prune_taxa(TopNOTUs, pvar.T0.ra)
pvar.T0.abund.melt <- psmelt(pvar.T0.abund)
pvar.T0.abund.melt$top.group <- "abund"
pvar.T0.abund.melt <- pvar.T0.abund.melt %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(pvar.T0.abund.melt, "../../top_abundance/output/pvar_T0_abund.csv")



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


#T0 only

pacu.T0 <- subset_samples(pacu, Time.Point == "T0")
pacu.T0.ra <- transform_sample_counts(pacu.T0, function(x) x/ sum(x))
TopNOTUs = names(sort(taxa_sums(pacu.T0.ra), TRUE)[1:10])
pacu.T0.abund = prune_taxa(TopNOTUs, pacu.T0.ra)
pacu.T0.abund.melt <- psmelt(pacu.T0.abund)
pacu.T0.abund.melt$top.group <- "abund"
pacu.T0.abund.melt <- pacu.T0.abund.melt %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(pacu.T0.abund.melt, "../../top_abundance/output/pacu_T0_abund.csv")


##Venn Diagram for T0 only
mcap.list <- unique(mcap.T0.abund.melt$otu_genus)
pcomp.list <- unique(pcomp.T0.abund.melt$otu_genus)
pvar.list <- unique(pvar.T0.abund.melt$otu_genus)
pacu.list <- unique(pacu.T0.abund.melt$otu_genus)

candidates <- list("Montipora capitata" = mcap.list, "Porites compressa" = pcomp.list,
                   "Pavona varians" = pvar.list, "Pocillopora acuta" = pacu.list)
library(VennDiagram)
venn <- venn.diagram(x = candidates, filename = NULL, fill = c("#E69F00", "#56B4E9", "#CC79A7","#009E73"), alpha = 0.5)
pdf(file = "../output/venn_abund_T0.pdf")
grid.draw(venn)
dev.off()
overlap <- calculate.overlap(candidates)
names(overlap) <- c("mcap pcomp pvar pacu", "mcap pcomp pvar", "mcap pcomp pacu", "mcap pvar pacu", "pcomp pvar pacu", "mcap pcomp", "mcap pvar ", 
                    "mcap pacu", "pcomp pvar", "pcomp  pacu", "pvar pacu", "mcap", "pcomp", "pvar", "pacu")
View(overlap)


##Can I do an endozoicomonas OTU one?

endo <- subset_taxa(bac.exp, Genus == "Endozoicomonas")
endo.melt <- psmelt(endo)
endo.melt <- filter(endo.melt, Time.Point == "T0")

mcap.endo <- filter(endo.melt, Species == "Montipora_capitata")
mcap.endo <- filter(mcap.endo, Abundance > 0)
mcap.list <- unique(mcap.endo$OTU)
pcomp.endo <- filter(endo.melt, Species == "Porites_compressa")
pcomp.endo <- filter(pcomp.endo, Abundance > 0)
pcomp.list <- unique(pcomp.endo$OTU)
pvar.endo <- filter(endo.melt, Species == "Pavona_varians")
pvar.endo <- filter(pvar.endo, Abundance > 0)
pvar.list <- unique(pvar.endo$OTU)
pacu.endo <- filter(endo.melt, Species == "Pocillopora_acuta")
pacu.endo <- filter(pacu.endo, Abundance > 0)
pacu.list <- unique(pacu.endo$OTU)

candidates <- list("Montipora capitata" = mcap.list, "Porites compressa" = pcomp.list,
                   "Pavona varians" = pvar.list, "Pocillopora acuta" = pacu.list)
library(VennDiagram)
venn <- venn.diagram(x = candidates, filename = NULL, fill = c("#E69F00", "#56B4E9", "#CC79A7","#009E73"), alpha = 0.5)
pdf(file = "../output/venn_endo_T0.pdf")
grid.draw(venn)
dev.off()
overlap <- calculate.overlap(candidates)
names(overlap) <- c("mcap pcomp pvar pacu", "mcap pcomp pvar", "mcap pcomp pacu", "mcap pvar pacu", "pcomp pvar pacu", "mcap pcomp", "mcap pvar ", 
                    "mcap pacu", "pcomp pvar", "pcomp  pacu", "pvar pacu", "mcap", "pcomp", "pvar", "pacu")
View(overlap)


