#### New R Script to streamline the 16S figure generation & analysis for ACE21 ####

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(plyr)

##NOTE: The working colours are as follows: 
# mcap "#E69F00" ; pcomp "#56B4E9" ; pvar "#CC79A7" ; pacu"#009E73"


#Check your working directory

getwd()

#Load 16S data in & extract the samples:

load("16s_phyloseq4HE.RData")
bac.exp <- subset_samples(Bac.seq, Type == "sample")


#Okay let's make some ordinations of the ambient time point

bac.ambient <- subset_samples(bac.exp, Time.Point == "T0")
bac.ambient.ra <- transform_sample_counts(bac.ambient, function(x) x/ sum(x)) 
bac.ambient.data <- as(sample_data(bac.ambient), "data.frame")

#Let's run a quick permanova to look at differences in bray curtis by species
bc.amb <- phyloseq::distance(bac.ambient.ra, method = "bray")
vegan::adonis2(bc.amb ~ Species, data = bac.ambient.data)
vegan::permutest(vegan::betadisper(bc.amb, bac.ambient.data$Species, type = "centroid"))
disp.amb <- vegan::betadisper(bc.amb, bac.ambient.data$Species, type = "centroid")


#1) 
#Extract distances from centroid to plot in ggplot
disp.amb.data <- as.data.frame(disp.amb$distances)
disp.amb.data$sample_name.1 <- rownames(disp.amb.data)
disp.amb.data <- merge(disp.amb.data, bac.ambient.data, by = "sample_name.1")

pdf(file = "bac_figures/output/organised_output/2_Distance_to_Centroid_ambient.pdf")
ggplot(disp.amb.data, aes(x = Species, y = `disp.amb$distances`)) +
  geom_boxplot(aes(color = Species)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  theme_classic() +
  ylab("Distance to Centroid") +
  xlab("\nHost Species") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 
dev.off()
  
##2)
#Build an ordination
ord.bc <- ordinate(bac.ambient.ra, "NMDS", "bray") 
vegan::stressplot(ord.bc)
scores(ord.bc, display="sites")
nmds.scores <- as.data.frame(scores(ord.bc)$sites)
nmds.scores$host <- bac.ambient.data$Species
nmds.scores$symbiont <- bac.ambient.data$Clade

hull_nmds <- nmds.scores %>%
  group_by(host) %>%
  slice(chull(NMDS1, NMDS2))

pdf(file = "bac_figures/output/organised_output/1_NMDS_bray_ambient.pdf")
ggplot(nmds.scores, aes(x = NMDS1, y= NMDS2)) + 
  geom_point(aes(colour = nmds.scores$host), size = 3) +
  geom_polygon(data = hull_nmds, aes(x=NMDS1,y=NMDS2,fill=host),alpha=0.30) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  labs(color = 'Host Species', fill = 'Host Species') +
  theme_classic() +
  #theme(legend.position = c(0.8, 0.75)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 
dev.off()


##Alpha diversity

alpha.amb <- phyloseq::estimate_richness(bac.ambient, measures = c("Observed", "Shannon"))
alpha.amb$Evenness <- alpha.amb$Shannon/log(alpha.amb$Observed)
alpha.amb$sample_name.1 <- rownames(alpha.amb)
alpha.amb <- merge(alpha.amb, bac.ambient.data, by = "sample_name.1")

pdf(file = "bac_figures/output/organised_output/3_Richness_ambient.pdf")
ggplot(alpha.amb, aes(x = Species, y = Observed)) +
  geom_boxplot(aes(color = Species)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab("Microbial Richness") +
  xlab("\nHost Species") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 
dev.off()


pdf(file = "bac_figures/output/organised_output/4_Evenness_ambient.pdf")
ggplot(alpha.amb, aes(x = Species, y = Evenness)) +
  geom_boxplot(aes(color = Species)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73")) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab("Microbial Community Evenness") +
  xlab("\nHost Species") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))
dev.off()


#Bar plots for phylum level and genus level 

#1 phylum level
percent.phy <- bac.ambient.ra %>% 
  tax_glom(taxrank = "Phylum") 
perc.melt.phy <- psmelt(percent.phy)

sum.phy <- ddply(perc.melt.phy, c("Phylum", "Species"), summarise,
                 N = length(Abundance),
                 mean = mean(Abundance),
                 sd = sd(Abundance), 
                 se = sd/sqrt(N)
)

ggplot(sum.phy, aes(x = Species, y = mean, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ggtitle("Communities by Phyla at Ambient Temperatures (T0)") +
  #guides(fill=NULL) +
  xlab("/nHost Species") +
  ylab("Mean Relative Abundance") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))


x = tapply(sum.phy$mean, sum.phy$Phylum, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x) #add top 10 to code below to set NAs

##Bacteria_Unclassified is not included as top 10 below bc we can't distinguish whether they all come from the same phylum
#Added # 11 = Acidobacteria, and moved unclassfied to Unassigned

sum.phy$Phylum = factor(as.character(sum.phy$Phylum), levels=names(x))
sum.phy$col_phylum <- sum.phy$Phylum
View(sum.phy$col_phylum)
#Set everything that is super low to NA so that we can call them "other"
sum.phy$col_phylum <- as.character(sum.phy$col_phylum)
sum.phy$col_phylum <- ifelse(is.na(sum.phy$col_phylum), 
                            'Unassigned', sum.phy$col_phylum)
sum.phy$col_phylum <- as.factor(sum.phy$col_phylum)

sum.phy$col_phylum[sum.phy$col_phylum != "Proteobacteria" &
                    sum.phy$col_phylum != "Bacteroidetes" & 
                    sum.phy$col_phylum != "Cyanobacteria" &
                    sum.phy$col_phylum != "Firmicutes" &
                    sum.phy$col_phylum != "Actinobacteria" &
                    sum.phy$col_phylum != "Epsilonbacteraeota" &
                    sum.phy$col_phylum != "Planctomycetes" &
                    sum.phy$col_phylum != "Spirochaetes" &
                    sum.phy$col_phylum != "Verrucomicrobia" &
                    sum.phy$col_phylum != "Acidobacteria" &
                    sum.phy$col_phylum != "Bacteria_Unassigned"] <- NA


levels(sum.phy$col_phylum)
# add new factor
sum.phy$col_phylum <- factor(sum.phy$col_phylum, levels = c(levels(sum.phy$col_phylum), "Other"))
# convert NAs to other
sum.phy$col_phylum[is.na(sum.phy$col_phylum)] = "Other"
sum.phy$col_phylum <- factor(x = sum.phy$col_phylum, levels = c("Other", "Proteobacteria", "Bacteroidetes",
                                                              "Cyanobacteria", "Firmicutes", "Actinobacteria",
                                                              "Epsilonbacteraeota", "Planctomycetes", "Spirochaetes", 
                                                              "Verrucomicrobia", "Acidobacteria"))
#Make a colour scheme
library(RColorBrewer)
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

#plot!
pdf(file = "bac_figures/output/organised_output/5_Phylum_bars_ambient.pdf")
ggplot(sum.phy, aes(x = Species, y = mean, fill = col_phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mycolors) +
  ylab("Relative Abundance") +
  xlab("\nHost Species") +
  ggtitle("Top Phyla at Ambient (T0)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#2
#How about with genus

percent.gen <- bac.ambient %>% 
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) 
perc.melt.gen <- psmelt(percent.gen)

sum.gen <- ddply(perc.melt.gen, c("Genus", "Species"), summarise,
                 N = length(Abundance),
                 mean = mean(Abundance),
                 sd = sd(Abundance), 
                 se = sd/sqrt(N)
)

ggplot(sum.gen, aes(x = Species, y = mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  ggtitle("Communities by Phyla at Ambient Temperatures (T0)") +
  guides(fill=FALSE) +
  xlab("\nHost Species") +
  ylab("Mean Relative Abundance") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))

##Let's just plot top 10 with all others set to "Other"
x = tapply(sum.gen$mean, sum.gen$Genus, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x) #add top 10 to code below to set NAs


sum.gen$Genus = factor(as.character(sum.gen$Genus), levels=names(x))
sum.gen$col_genus<- sum.gen$Genus
#Set everything that is super low to NA so that we can call them "other"
sum.gen$col_genus <- as.character(sum.gen$col_genus)
sum.gen$col_genus <- ifelse(is.na(sum.gen$col_genus), 
                             'Unassigned', sum.gen$col_genus)
sum.gen$col_genus<- as.factor(sum.gen$col_genus)

sum.gen$col_genus[sum.gen$col_genus!= "Endozoicomonas" &
                    sum.gen$col_genus != "P3OB-42_ge" & 
                    sum.gen$col_genus != "Candidatus_Amoebophilus" &
                    sum.gen$col_genus != "MBIC10086" &
                    sum.gen$col_genus != "Fulvivirga" &
                    sum.gen$col_genus != "Mastigocoleus_BC008" &
                    sum.gen$col_genus != "Thalassobius" &
                    sum.gen$col_genus != "Pelagibius" &
                    sum.gen$col_genus != "Alteromonas" &
                    sum.gen$col_genus != "Alkalispirochaeta" ] <- NA # &
                    #sum.gen$col_genus != "Unassigned"] <- NA


levels(sum.gen$col_genus)
# add new factor
sum.gen$col_genus <- factor(sum.gen$col_genus, levels = c(levels(sum.gen$col_genus), "Other"))
# convert NAs to other
sum.gen$col_genus[is.na(sum.gen$col_genus)] = "Other"
sum.gen$col_genus <- factor(x = sum.gen$col_genus, levels = c("Other", "Endozoicomonas", "P3OB-42_ge","Candidatus_Amoebophilus",
                                                              "MBIC10086", "Fulvivirga", "Mastigocoleus_BC008",
                                                              "Thalassobius", "Pelagibius", "Alteromonas", 
                                                              "Alkalispirochaeta"))
#Make a colour scheme
library(RColorBrewer)
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

#plot!
pdf(file = "bac_figures/output/organised_output/6_Genus_bars_ambient.pdf")
ggplot(sum.gen, aes(x = Species, y = mean, fill = col_genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mycolors) +
  ylab("Relative Abundance") +
  xlab("\nHost Species") +
  ggtitle("Top Genera at Ambient (T0)")+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


