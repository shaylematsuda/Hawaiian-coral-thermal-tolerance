##Creating the Figure of Top Abundant + Top Core + Top Differential by Host species
##For Matsuda et al. 2021

#load libraries
library(phyloseq)
library(ggplot2)

#Set working directory
setwd("~/Documents/OSUDocs/Projects/ACE21/Hawaiian-coral-thermal-tolerance/Data/16s_analysis/")

#bring in data from ANCOM
mcap.ancom.melt <- read.csv("ancom/output/mcap.ancom.csv")
pcomp.ancom.melt <- read.csv("ancom/output/pcomp_ancom.csv")
pvar.ancom.melt <- read.csv("") #So there are no significantly differentially abundant taxa for p var by treatment
pacu.ancom.melt <- read.csv("")

#bring in data from CORE
mcap.core.melt <- read.csv("core_microbiome/output/mcap_core.csv")
pcomp.core.melt <- read.csv("core_microbiome/output/pcomp_core.csv")
pvar.core.melt <- read.csv("core_microbiome/output/pvar_core.csv")
pacu.core.melt <- read.csv("core_microbiome/output/pacu_core.csv")

#bring in the data from top abundant taxa
mcap.abund.melt <- read.csv("top_abundance/output/mcap_abund.csv")
pcomp.abund.melt <-read.csv("top_abundance/output/pcomp_abund.csv")
pvar.abund.melt <- read.csv("top_abundance/output/pvar_abund.csv")
pacu.abund.melt <- read.csv("top_abundance/output/pacu_abund.csv")


##1 Montipora Capitata
#First this is a work around to get a list of the samples in the right order for plotting 
#in order of: time, treatment, clade
mcap.toget <- mcap.core.melt[order(mcap.core.melt$Time.Point, mcap.core.melt$Treatment, mcap.core.melt$Clade),]
mcap.toget <- filter(mcap.toget, OTU == "Otu00001")
mcap.list <- as.list(mcap.toget$Sample)

#merge all data
mcap.top <- rbind(mcap.ancom.melt, mcap.core.melt, mcap.abund.melt)

p.mcap <- ggplot(mcap.top, aes(x = Sample, y = Genus, size = Abundance, color = Time.Point, shape = Treatment)) +
  geom_point() +
  scale_shape_manual(values = c(16,18)) +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("sandybrown", "tomato4","darkslategray4")) +
  facet_grid(top.group ~ .) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p.mcap$data$Sample <- as.factor(p.mcap$data$Sample)
p.mcap$data$Sample <- factor(p.mcap$data$Sample, levels = mcap.list)
print(p.mcap)
ggsave("bac_figures/output/mcap_top_bubble.pdf")
##I also saved this as png from plot zoom because they are easier to read 

##2 Porites compressa
#First this is a work around to get a list of the samples in the right order for plotting 
#in order of: time, treatment, clade
pcomp.toget <- pcomp.ancom.melt[order(pcomp.ancom.melt$Time.Point, pcomp.ancom.melt$Treatment, pcomp.ancom.melt$Clade),]
pcomp.toget <- filter(pcomp.toget, OTU == "Otu00073")
pcomp.list <- as.list(pcomp.toget$Sample)

#merge all data
pcomp.top <- rbind(pcomp.ancom.melt, pcomp.core.melt, pcomp.abund.melt)

p.pcomp <- ggplot(pcomp.top, aes(x = Sample, y = OTU, size = Abundance, color = Time.Point, shape = Treatment)) +
  geom_point() +
  scale_shape_manual(values = c(16,18)) +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("sandybrown", "tomato4","darkslategray4")) +
  facet_grid(top.group ~ .) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p.pcomp$data$Sample <- as.factor(p.pcomp$data$Sample)
p.pcomp$data$Sample <- factor(p.pcomp$data$Sample, levels = pcomp.list)
print(p.pcomp)
ggsave("bac_figures/output/pcomp_top_bubble.pdf")


##3 Pavona varians
#First this is a work around to get a list of the samples in the right order for plotting 
#in order of: time, treatment, clade
pvar.toget <- pvar.core.melt[order(pvar.core.melt$Time.Point, pvar.core.melt$Treatment, pvar.core.melt$Clade),]
pvar.toget <- filter(pvar.toget, OTU == "Otu00001")
pvar.list <- as.list(pvar.toget$Sample)

#merge all data
pvar.top <- rbind(pvar.core.melt, pvar.abund.melt)

p.pvar<- ggplot(pvar.top, aes(x = Sample, y = OTU, size = Abundance, color = Time.Point, shape = Treatment)) +
  geom_point() +
  scale_shape_manual(values = c(16,18)) +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("sandybrown", "tomato4","darkslategray4")) +
  facet_grid(top.group ~ .) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p.pvar$data$Sample <- as.factor(p.pvar$data$Sample)
p.pvar$data$Sample <- factor(p.pvar$data$Sample, levels = pvar.list)
print(p.pvar)
ggsave("bac_figures/output/pvar_top_bubble.pdf")



##4 Pocillopora acuta
#First this is a work around to get a list of the samples in the right order for plotting 
#in order of: time, treatment, clade
pacu.core.melt
levels(pacu.core.melt$Time.Point)[levels(pacu.core.melt$Time.Point)=="F1"] <- "T1"
pacu.toget <- pacu.core.melt[order(pacu.core.melt$Time.Point, pacu.core.melt$Treatment, pacu.core.melt$Clade),]
pacu.toget <- filter(pacu.toget, OTU == "Otu00019")
pacu.list <- as.list(pacu.toget$Sample)

#merge all data
pacu.top <- rbind(pacu.ancom.melt, pacu.core.melt, pacu.abund.melt)

p.pacu <- ggplot(pacu.top, aes(x = Sample, y = OTU, size = Abundance, color = Time.Point, shape = Treatment)) +
  geom_point() +
  scale_shape_manual(values = c(16,18)) +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("sandybrown", "tomato4","darkslategray4")) +
  facet_grid(top.group ~ .) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p.pacu$data$Sample <- as.factor(p.pacu$data$Sample)
p.pacu$data$Sample <- factor(p.pacu$data$Sample, levels = pacu.list)
print(p.pacu)
ggsave("bac_figures/output/pacu_top_bubble.pdf")

