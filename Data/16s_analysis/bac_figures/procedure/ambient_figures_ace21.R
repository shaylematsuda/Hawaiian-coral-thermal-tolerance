#### New R Script to streamline the 16S figure generation & analysis for ACE21 ####

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(plyr)
library(vegan)

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
#Df SumOfSqs      R2      F Pr(>F)    
#Species   3   5.6919 0.32431 6.0797  0.001 ***
#  Residual 38  11.8586 0.67569                  
#Total    41  17.5505 1.00000   
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
#Re-order the coral species so that it is M cap, P comp, P var, P acu
alpha.amb$Species <- factor(alpha.amb$Species,
                       levels = c('Montipora_capitata','Porites_compressa', "Pavona_varians", "Pocillopora_acuta"),ordered = TRUE)
#This means the order of colors should be: "#E69F00", "#009E73" ,"#56B4E9","#CC79A7"

pdf(file = "bac_figures/output/organised_output/3_Richness_ambient.pdf")
ggplot(alpha.amb, aes(x = Species, y = Observed)) +
  geom_boxplot(aes(fill = Species)) +
  scale_fill_manual(values=c("#E69F00", "#009E73" ,"#56B4E9","#CC79A7")) +
  scale_colour_manual(values=c("#E69F00", "#009E73" ,"#56B4E9","#CC79A7")) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab("Microbial Richness") +
  xlab("\nHost Species") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 
dev.off()


pdf(file = "bac_figures/output/organised_output/4_Evenness_ambient.pdf")
ggplot(alpha.amb, aes(x = Species, y = Evenness)) +
  geom_boxplot(aes(fill = Species)) +
  scale_fill_manual(values=c("#E69F00", "#009E73" ,"#56B4E9","#CC79A7")) +
  scale_colour_manual(values=c("#E69F00", "#009E73" ,"#56B4E9","#CC79A7")) +
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
#Might be an issue with saying "Species" as host species

#This code works: 

colnames(sample_data(bac.ambient.ra))[colnames(sample_data(bac.ambient.ra)) == "Species"] <- "host_species"
head(sample_data(bac.ambient.ra))

percent.gen <- bac.ambient.ra %>% 
  tax_glom(taxrank = "Genus") %>%
  merge_samples("host_species") %>%  
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()

percent.gen$Genus<- as.character(percent.gen$Genus) #convert to character
percent.gen$Genus[percent.gen$Abundance < 0.02] <- "< 2% abund."

#Apparently changes the merge factor to "Sample" column - so we need to make sure these are in the right order: 
percent.gen$Sample <- as.factor(percent.gen$Sample)
percent.gen$Sample <- factor(percent.gen$Sample, levels = c("Pocillopora_acuta","Pavona_varians",  "Porites_compressa", "Montipora_capitata"))


p.gen.amb <- ggplot(percent.gen, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  ggtitle("Communities by Genus at Ambient Temperatures (T0)") +
  guides() +
  scale_fill_manual(values = c('#AAAAAB','#2B8CBC',  '#9ABCAD','#D6BAD9', '#F6BD91', '#9FD7F5', '#AA76B3', 
                               '#70BADD',  '#AA7552','#C97D67',  '#8D6373', '#DBF1F9', 
                               '#476269')) +
  xlab("\nHost Species") +
  scale_x_discrete(labels= c("P.acuta", "P. varians", "P. compressa", "M. capitata")) +
  ylab("Mean Relative Abundance") +
  coord_flip() +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))


ggsave("bac_figures/output/organised_output/6b_Genus_bars_ambient_coord.pdf", plot = p.gen.amb)


#Not quite built the same way as the phylum level plot

#How about building a Phylum plot the same way then? 

percent.phy <- bac.ambient.ra %>% 
  tax_glom(taxrank = "Phylum") %>%
  merge_samples("host_species") %>%  
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()

percent.phy$Phylum<- as.character(percent.phy$Phylum) #convert to character
percent.phy$Phylum[percent.phy$Abundance < 0.01] <- "< 1% abund."

percent.phy$Sample <- as.factor(percent.phy$Sample)
percent.phy$Sample <- factor(percent.phy$Sample, levels = c("Pocillopora_acuta","Pavona_varians",  "Porites_compressa", "Montipora_capitata"))

#species_names <- c("M. capitata", "P. varians", "P. acuta", "P. compressa")

p.phy.amb <- ggplot(percent.phy, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  ggtitle("Communities by Phylum at Ambient Temperatures (T0)") +
  guides() +
  scale_fill_manual(values = c('#808080','#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032',
                               '#C2B280','#0067A5',  '#E68FAC', '#654522')) +
  xlab("\nHost Species") +
  scale_x_discrete(labels= c("P. acuta","P. varians", "P. compressa", "M. capitata")) +
  ylab("Mean Relative Abundance") +
  coord_flip() +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))


ggsave("bac_figures/output/organised_output/5b_Phylum_bars_ambient_coord.pdf", plot = p.phy.amb)



##Can we put them on the same graph on top of each other? 
install.packages("patchwork")
library(patchwork)
 # Combine the plots
gen.phy.wrap <- wrap_plots(p.phy.amb, p.gen.amb, ncol = 1, nrow = 2)

ggsave("bac_figures/output/organised_output/5_6_bars_ambient_coord_wrap.pdf", plot = gen.phy.wrap)



# Genus level plot currently deprecated version below here
#Commented out, but could be useful in future:


#percent.gen <- bac.ambient.ra %>% 
#  tax_glom(taxrank = "Genus") %>% 
#  transform_sample_counts(function(x) {x/sum(x)}) %>%
#  psmelt()
#View(percent.gen)
#View(tax_table(bac.ambient))
#perc.melt.gen <- psmelt(percent.gen)


#ggplot(percent.gen, aes(x = Species, y = Abundance, fill = Genus))+
#  geom_bar(stat = "identity") +
#  ggtitle("Communities by Genus at Ambient Temperatures (T0)") +
#  guides(fill="none") +
#  xlab("\nHost Species") +
#  ylab("Mean Relative Abundance") +
#  theme_classic() +
#  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))


#sum.gen <- ddply(percent.gen, c("Genus", "Species"), summarise,
#                 N = length(Abundance),
#                 mean = mean(Abundance),
#                 sd = sd(Abundance), 
#                 se = sd/sqrt(N)
#)

#ggplot(sum.gen, aes(x = Species, y = mean, fill = Genus)) +
#  geom_bar(stat = "identity") +
#  ggtitle("Communities by Phyla at Ambient Temperatures (T0)") +
#  guides(fill="none") +
#  xlab("\nHost Species") +
#  ylab("Mean Relative Abundance") +
#  theme_classic() +
#  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))

##Let's just plot top 10 with all others set to "Other"
#x = tapply(sum.gen$mean, sum.gen$Genus, function(x) max(x))
#x = sort(x, TRUE) #sort phyla 
#View(x) #add top 10 to code below to set NAs


#sum.gen$Genus = factor(as.character(sum.gen$Genus), levels=names(x))
#sum.gen$col_genus<- sum.gen$Genus
#Set everything that is super low to NA so that we can call them "other"
#sum.gen$col_genus <- as.character(sum.gen$col_genus)
#sum.gen$col_genus <- ifelse(is.na(sum.gen$col_genus), 
#                             'Unassigned', sum.gen$col_genus)
#sum.gen$col_genus<- as.factor(sum.gen$col_genus)

#sum.gen$col_genus[sum.gen$col_genus!= "Endozoicomonas" &
#                    sum.gen$col_genus != "P3OB-42_ge" & 
#                    sum.gen$col_genus != "Candidatus_Amoebophilus" &
#                    sum.gen$col_genus != "MBIC10086" &
#                    sum.gen$col_genus != "Fulvivirga" &
#                    sum.gen$col_genus != "Mastigocoleus_BC008" &
#                    sum.gen$col_genus != "Thalassobius" &
#                    sum.gen$col_genus != "Pelagibius" &
#                    sum.gen$col_genus != "Alteromonas" &
#                    sum.gen$col_genus != "Alkalispirochaeta" ] <- NA # &
                    #sum.gen$col_genus != "Unassigned"] <- NA


#levels(sum.gen$col_genus)
# add new factor
#sum.gen$col_genus <- factor(sum.gen$col_genus, levels = c(levels(sum.gen$col_genus), "Other"))
# convert NAs to other
#sum.gen$col_genus[is.na(sum.gen$col_genus)] = "Other"
#sum.gen$col_genus <- factor(x = sum.gen$col_genus, levels = c("Other", "Endozoicomonas", "P3OB-42_ge","Candidatus_Amoebophilus",
#                                                              "MBIC10086", "Fulvivirga", "Mastigocoleus_BC008",
#                                                              "Thalassobius", "Pelagibius", "Alteromonas", 
#                                                              "Alkalispirochaeta"))
#Make a colour scheme
#library(RColorBrewer)
#nb.cols <- 11
#mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

#plot!
#pdf(file = "bac_figures/output/organised_output/6_Genus_bars_ambient.pdf")
#ggplot(sum.gen, aes(x = Species, y = mean, fill = col_genus)) +
#  geom_bar(stat = "identity") +
#  scale_fill_manual(values=mycolors) +
#  ylab("Relative Abundance") +
#  xlab("\nHost Species") +
#  ggtitle("Top Genera at Ambient (T0)")+
#  theme_classic() +
#  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#dev.off()


