###Code for core microbiome analyses for Matsuda et al. 2021###
##Purpose: to look for shared taxa among species

#Load libraries
library(phyloseq)
library(plyr)
library(microbiome)

#Set working directory & bring in RData
#source the ANCOM functions R script
setwd("~/Documents/OSUDocs/Projects/ACE21/Hawaiian-coral-thermal-tolerance/Data/16s_analysis/")
load("16s_phyloseq4HE.RData")

#subset to just experimental samples
bac.exp <- subset_samples(Bac.seq, Type == "sample")

##Check out the core for all samples at T0 & renove OTUs with 0 abundance
t0 <- subset_samples(bac.exp, Time.Point == "T0")
t0 <- prune_taxa(taxa_sums(t0) > 0, t0)

##Core Analysis
t0.core <- core(t0, detection = 0.0001, prevalence = .78)
t0.core.taxa <- tax_table(t0.core)
View(as.data.frame(t0.core.taxa))
##Only OTU00001 Endozoicomonas is a core within 78% of samples. Nothing above this level

#Trial with a tax glom at family??
t0.fam <- tax_glom(t0, "Family", NArm = TRUE)
t0.fam.core <- core(t0.fam, detection = 0.0001, prevalence = .88)
t0.fam.core.taxa <- tax_table(t0.fam.core)
View(as.data.frame(t0.fam.core.taxa))
##OTU00693 Flavibacteraceae
#OTU00073 Rhodobacteraceae



##Okay what about within each species?

#First up M. capitata
##Here we are using a prevlance of 60% to get 5 "core" taxa - top most shared 
mcap <- subset_samples(bac.exp, Species == "Montipora_capitata")
mcap.ra <- transform_sample_counts(mcap, function(x) x/ sum(x))
mcap.core <- core(mcap.ra, detection = 0.0001, prevalence = 0.6)
mcap.core.taxa <- tax_table(mcap.core)
View(as.data.frame(mcap.core.taxa))
mcap.core.melt <- psmelt(mcap.core)
mcap.core.melt$top.group <- "core"
write.csv(mcap.core.melt, "core_microbiome/output/mcap_core.csv")

## There are 5

##Plotting
#check # of OTUs
mcap.core.melt$OTU <- as.factor(mcap.core.melt$OTU)
levels(mcap.core.melt$OTU) 
#re-order according to Time point and clade
mcap.core.ordered <- mcap.core.melt[order(mcap.core.melt$Time.Point, mcap.core.melt$Clade),]


p.mcap <- ggplot(mcap.core.ordered, aes(x = Sample, y = OTU, size = Abundance, color = Treatment)) +
  geom_point() +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("dodgerblue4", "darkorange2","indianred3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p.mcap$data$Sample <- as.factor(p.mcap$data$Sample)
p.mcap$data$Sample <- factor(p.mcap$data$Sample, levels = c("SN38", "SN21", "SN9", "SN94", "S1753b",
                                                            "SN82", "S1753","SN33", "SN64",
                                                            "SN72", "S419", "S848", "SN192", "SN184", 
                                                            "SN49", "S153", "S1791", "S3479", "S3472",
                                                            "S3482", "SN92", "SN45", "SN172", "S797", "SN190", "S3469", "S179",
                                                            "S1133", "SN150", "S1226", "S377","S3473",  "S335", 
                                                            "S1437", "S3422", "S3204", "SN42", "SN152", "SN110", "SN96", 
                                                            "S1104","S3470", "S2284", "S619"))
print(p.mcap)
ggsave("core_microbiome/output/mcap_core_bubble.pdf")


##Porites compressa
##Use 50% to get at least 5 top shared taxa
pcomp <- subset_samples(bac.exp, Species == "Porites_compressa")
pcomp.ra <- transform_sample_counts(pcomp, function(x) x/ sum(x))
pcomp.core <- core(pcomp.ra, detection = 0.0001, prevalence = 0.5)
pcomp.core.taxa <- tax_table(pcomp.core)
View(as.data.frame(pcomp.core.taxa))
pcomp.core.melt <- psmelt(pcomp.core)
pcomp.core.melt$top.group <- "core"
write.csv(pcomp.core.melt, "core_microbiome/output/pcomp_core.csv")

##Plotting
#check # of OTUs
pcomp.core.melt$OTU <- as.factor(pcomp.core.melt$OTU)
levels(pcomp.core.melt$OTU) 
#re-order according to Time point and clade
pcomp.core.ordered <- pcomp.core.melt[order(pcomp.core.melt$Time.Point, pcomp.core.melt$Clade),]


p.pcomp <- ggplot(pcomp.core.ordered, aes(x = Sample, y = OTU, size = Abundance, color = Time.Point)) +
  geom_point() +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("dodgerblue4", "darkorange2","indianred3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p.pcomp$data$Sample <- as.factor(p.pcomp$data$Sample)
p.pcomp$data$Sample <- factor(p.pcomp$data$Sample, levels = c("SN3", "SN19", "S74", "S31", "SN7", "SN107", "SN50",
                                                              "SN47", "SN25", "SN70b", "S650", "SN3340", "S1110", "SN117", "S199",
                                                              "S284", "SN149", "S3458", "SN37", "S3514", "SN157", "S3340", "SN134",
                                                              "S117", "S226", "SN78", "S634", "SN124", "S3488",
                                                              "SN181", "S2286", "S946", "S809", "SN161", "S3475", "SN115", "SN122", 
                                                              "S3476b", "SN105", "S3476", "S3491", "SN77", "S1502","S3485", "S3517",
                                                              "S3480", "SN205", "SN98", "S3337", "S579", "S3467", "S3442", "SN131", 
                                                              "S972", "S3444", "SN211", "SN189", "SN55"))
print(p.pcomp)
ggsave("core_microbiome/output/pcomp_core_bubble.pdf")


##Pavona varians
##Use 45% to get at 4 top shared taxa
pvar <- subset_samples(bac.exp, Species == "Pavona_varians")
pvar.ra <- transform_sample_counts(pvar, function(x) x/ sum(x))
pvar.core <- core(pvar.ra, detection = 0.0001, prevalence = 0.45)
pvar.core.taxa <- tax_table(pvar.core)
View(as.data.frame(pvar.core.taxa))
pvar.core.melt <- psmelt(pvar.core)
pvar.core.melt$top.group <- "core"
write.csv(pvar.core.melt, "core_microbiome/output/pvar_core.csv")

##Plotting
#check # of OTUs
pvar.core.melt$OTU <- as.factor(pvar.core.melt$OTU)
levels(pvar.core.melt$OTU) 
#re-order according to Time point and clade
pvar.core.ordered <- pvar.core.melt[order(pvar.core.melt$Time.Point, pvar.core.melt$Clade),]
View(pvar.core.ordered)

p.pvar <- ggplot(pvar.core.ordered, aes(x = Sample, y = OTU, size = Abundance, color = Time.Point)) +
  geom_point() +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("dodgerblue4", "darkorange2","indianred3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p.pvar$data$Sample <- as.factor(p.pvar$data$Sample)
p.pvar$data$Sample <- factor(p.pvar$data$Sample, levels = c("S914", "S488", "S148","SN35", "SN35b", "S263", "S616",
                                                             "SN84", "S231", "S208", "S208b", "S391", "S280b", "S1982", 
                                                             "S1108", "S3478", "S3246", "SN108b", "S1068", "S664", 
                                                             "S1943", "S1144", "S1441", "SN222", "S761", "S3137",
                                                             "S249", "S2629", "S2911", "S436", "S1254", "S490", "S2001", 
                                                             "S554", "S1857", "S2849", "S2731", "SN217", "S1803", "S3383", 
                                                             "S955", "S775", "S1037", "S1526", "S471", "SN191",
                                                            "S2663", "S2214", "S1413", "S1688", "S422","S2328", "S2136", 
                                                            "S832", "S1556", "S2948", "S420", "S1166"))
print(p.pvar)
ggsave("core_microbiome/output/pvar_core_bubble.pdf")



##Pocillopora acuta
##Use 57% to get at 4 top shared taxa
pacu <- subset_samples(bac.exp, Species == "Pocillopora_acuta")
pacu.ra <- transform_sample_counts(pacu, function(x) x/ sum(x))
pacu.core <- core(pacu.ra, detection = 0.0001, prevalence = 0.57)
pacu.core.taxa <- tax_table(pacu.core)
View(as.data.frame(pacu.core.taxa))
pacu.core.melt <- psmelt(pacu.core)
pacu.core.melt$top.group <- "core"
write.csv(pacu.core.melt, "core_microbiome/output/pacu_core.csv")

##Plotting
#check # of OTUs
pacu.core.melt$OTU <- as.factor(pacu.core.melt$OTU)
levels(pacu.core.melt$OTU) 
#re-order according to Time point and clade
pacu.core.ordered <- pacu.core.melt[order(pacu.core.melt$Time.Point, pacu.core.melt$Clade),]
#Just use on OTU to get the correct order of sampleIDs and make into a list
pacu.core.toget <- filter(pacu.core.ordered, OTU == "Otu00001")
View(pacu.core.toget)
pacu.list <- as.list(pacu.core.toget$Sample)

p.pacu <- ggplot(pacu.core.ordered, aes(x = Sample, y = OTU, size = Abundance, color = Time.Point)) +
  geom_point() +
  scale_size_area(max_size = 6, breaks = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )) +
  coord_fixed(ratio = 0.9) +
  scale_color_manual(values=c("dodgerblue4", "darkorange2","indianred3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p.pacu$data$Sample <- as.factor(p.pacu$data$Sample)
p.pacu$data$Sample <- factor(p.pacu$data$Sample, levels = list)
print(p.pacu)
ggsave("core_microbiome/output/pacu_core_bubble.pdf")


