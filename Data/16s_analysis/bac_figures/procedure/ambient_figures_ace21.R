#### New R Script to streamline the 16S figure generation & analysis for ACE21 ####

#Load libraries
check

library(phyloseq)
library(tidyverse)
##NOTE: The working colours are as follows: 
# mcap "#E69F00" ; pcomp "#56B4E9" ; pvar "#CC79A7" ; pacu"#009E73"


#Check your working directory

getwd()
#want to be in /Data/16s_analysis/bac_figures/procedure/
#Load 16S data in & extract the samples:

load("../../16s_phyloseq4HE.RData")
bac.exp <- subset_samples(Bac.seq, Type == "sample")

#We have an issue with rank names
rank_names(bac.exp)

#Let's rename them: 
bac.exp.1 <- bac.exp %>% speedyseq::rename_tax_table(Kingdom = Domain) %>%
  speedyseq::rename_tax_table(Order2 = Class) %>%
  speedyseq::rename_tax_table(Class = Order) %>%
  speedyseq::rename_tax_table(Order = Order2)
rank_names(bac.exp.1)

#Okay let's make some ordinations of the ambient time point

bac.ambient <- subset_samples(bac.exp.1, Time.Point == "T0")
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
disp.amb.data$Species <- factor(disp.amb.data$Species,
                           levels = c('Montipora_capitata','Porites_compressa', "Pavona_varians", "Pocillopora_acuta"),ordered = TRUE)




pdf(file = "../output/organised_output/2_Distance_to_Centroid_ambient.pdf")
ggplot(disp.amb.data, aes(x = Species, y = `disp.amb$distances`)) +
  geom_boxplot(aes(color = Species)) +
  scale_fill_manual(values=c("#E69F00","#56B4E9","#CC79A7", "#009E73")) +
  scale_colour_manual(values=c("#E69F00","#56B4E9","#CC79A7", "#009E73")) +
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

nmds.scores$host <- factor(nmds.scores$host,
                                levels = c('Montipora_capitata','Porites_compressa', "Pavona_varians", "Pocillopora_acuta"),ordered = TRUE)



hull_nmds <- nmds.scores %>%
  group_by(host) %>%
  slice(chull(NMDS1, NMDS2))


pdf(file = "../output/organised_output/1_NMDS_bray_ambient.pdf")
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

pdf(file = "../output/organised_output/3_Richness_ambient.pdf")
ggplot(alpha.amb, aes(x = Species, y = Observed)) +
  geom_boxplot(aes(fill = Species)) +
  scale_fill_manual(values=c("#E69F00","#56B4E9","#CC79A7", "#009E73")) +
  scale_colour_manual(values=c("#E69F00","#56B4E9","#CC79A7", "#009E73")) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylab("Microbial Richness") +
  xlab("\nHost Species") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 
dev.off()


pdf(file = "../output/organised_output/4_Evenness_ambient.pdf")
ggplot(alpha.amb, aes(x = Species, y = Evenness)) +
  geom_boxplot(aes(fill = Species)) +
  scale_fill_manual(values=c("#E69F00","#56B4E9","#CC79A7", "#009E73" )) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9","#CC79A7", "#009E73")) +
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

sum.gen <- plyr::ddply(percent.gen, c("Genus", "host_species"), summarise,
                 N = length(Abundance),
                 mean = mean(Abundance),
                 sd = sd(Abundance), 
                 se = sd/sqrt(N)
)
View(sum.gen)

percent.gen$Genus<- as.character(percent.gen$Genus) #convert to character
percent.gen$Genus[percent.gen$Abundance < 0.02] <- "< 2% abund."

#Apparently changes the merge factor to "Sample" column - so we need to make sure these are in the right order: 
percent.gen$Sample <- as.factor(percent.gen$Sample)
percent.gen$Sample <- factor(percent.gen$Sample, levels = c("Montipora_capitata", "Porites_compressa", "Pavona_varians", "Pocillopora_acuta"))


p.gen.amb <- ggplot(percent.gen, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  #ggtitle("Communities by Genus at Ambient Temperatures (T0)") +
  guides() +
  scale_fill_manual(values = c('#AAAAAB','#2B8CBC',  '#9ABCAD','#D6BAD9', '#F6BD91', '#9FD7F5', '#AA76B3', 
                               '#70BADD',  '#AA7552','#C97D67',  '#8D6373', '#DBF1F9', 
                               '#476269'),
                    labels = c("<2% abundance", "Alphaproteobacteria (unclassified)", "Bacteria (unclassified)", 
                              expression(italic("Ca.")*" Amoebophilus"), "Cytophagales (unclassified)", expression(italic("Endozoicomonas")), 
                              expression(italic("Fulvivirga")), "Gammaproteobacteria (unclassified)", expression(italic("Mastigocoleus")*" BC008"),
                              "Cyanobacterium MBIC10086", "Myxococcales P3OB-42_ge", "Rhodobacteriaceae (unclassified)",
                              "Thalassobaculales (uncultured)")) +
  #xlab("\nHost Species") +
  xlab("") +
  scale_x_discrete(labels= c("M. capitata", "P. compressa", "P. varians", "P.acuta")) +
  ylab("Mean Relative Abundance") +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  #coord_flip() + #this flips the plot on its side
  theme_classic() +
  #theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.line.x = element_line(linewidth = 0.15), axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8, face = 'bold.italic')) +
  theme(axis.line.y = element_line(linewidth = 0.15), axis.title.y = element_text(size = 12), axis.ticks.y = element_blank()) +
  theme( legend.box.margin = margin(t = 50), legend.text = element_text(size = 8))


ggsave("../output/organised_output/6b_Genus_bars_ambient_coord.pdf", plot = p.gen.amb)


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
percent.phy$Sample <- factor(percent.phy$Sample, levels = c("Montipora_capitata", "Porites_compressa", "Pavona_varians", "Pocillopora_acuta"))

#species_names <- c("M. capitata", "P. varians", "P. acuta", "P. compressa")

p.phy.amb <- ggplot(percent.phy, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  #ggtitle("Communities by Phylum at Ambient Temperatures (T0)") +
  guides() +
  scale_fill_manual(values = c('#808080','#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032',
                               '#C2B280','#0067A5',  '#E68FAC', '#654522'),  
                    labels = c("<1% abundance", "Actinobacteriota", "Bacteria (unclassified)", "Bacteroidota",
                               "Cyanobacteriota", "Epsilonbacteraeota", "Bacillota", "Planctomycetota", 
                               "Pseudomonadota", "Spirochaetota", "Verrucomicrobiota")) +
  #xlab("\nHost Species") +
  xlab("") +
  scale_x_discrete(labels= c("M. capitata", "P. compressa", "P. varians", "P. acuta")) +
  ylab("Mean Relative Abundance") +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  #coord_flip() + #this flips the plot on its side
  theme_classic() +
  #theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.line.x = element_line(linewidth = 0.15) , axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8, face = 'bold.italic')) +
  theme(axis.line.y = element_line(linewidth = 0.15), axis.title.y = element_text(size = 12), axis.ticks.y = element_blank()) +
  theme(legend.box.margin = margin(t = 20), legend.text = element_text(size = 8))
  

ggsave("../output/organised_output/5b_Phylum_bars_ambient_coord.pdf", plot = p.phy.amb)

#3 Can we make an order level plot??
percent.order <- bac.ambient.ra %>% 
  tax_glom(taxrank = "Order") 
perc.melt.order <- psmelt(percent.order)

sum.order <- plyr::ddply(perc.melt.order, c("Order", "Species"), summarise,
                 N = length(Abundance),
                 mean = mean(Abundance),
                 sd = sd(Abundance), 
                 se = sd/sqrt(N)
)

ggplot(sum.order, aes(x = Species, y = mean, fill = Order)) +
  geom_bar(stat = "identity") +
  ggtitle("Communities by Order at Ambient Temperatures (T0)") +
  #guides(fill=NULL) +
  xlab("/nHost Species") +
  ylab("Mean Relative Abundance") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))


x = tapply(sum.order$mean, sum.order$Order, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x) #add top 10 to code below to set NAs

##Bacteria_Unclassified is not included as top 10 below bc we can't distinguish whether they all come from the same phylum
#Added # 11 = Rhizobiales, and moved unclassified to Unassigned

sum.order$Order = factor(as.character(sum.order$Order), levels=names(x))
sum.order$col_order <- sum.order$Order
View(sum.order$col_order)
#Set everything that is super low to NA so that we can call them "other"
sum.order$col_order <- as.character(sum.order$col_order)
sum.order$col_order <- ifelse(is.na(sum.order$col_order), 
                             'Unassigned', sum.order$col_order)
sum.order$col_order <- as.factor(sum.order$col_order)

sum.order$col_order[sum.order$col_order != "Oceanospirillales" &
                     sum.order$col_order != "Myxococcales" & 
                     sum.order$col_order != "Cytophagales" &
                     sum.order$col_order != "Phormidesmiales" &
                     sum.order$col_order != "Rhodobacterales" &
                     sum.order$col_order != "Flavobacteriales" &
                     sum.order$col_order != "Alteromonadales" &
                     sum.order$col_order != "Rickettsiales" &
                     sum.order$col_order != "Cellvibrionales" &
                     sum.order$col_order != "Rhizobiales"] <- NA


levels(sum.order$col_order)
# add new factor
sum.order$col_order <- factor(sum.order$col_order, levels = c(levels(sum.order$col_order), "Other"))
# convert NAs to other
sum.order$col_order[is.na(sum.order$col_order)] = "Other"
sum.order$col_order <- factor(x = sum.order$col_order, levels = c("Other", "Rhizobiales",
                                                                  "Cellvibrionales", "Rickettsiales", "Alteromonadales",
                                                                  "Flavobacteriales", "Rhodobacterales" , "Phormidesmiales", 
                                                                  "Cytophagales", "Myxococcales","Oceanospirillales"))
#Make a colour scheme
library(RColorBrewer)
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

#plot!
ggplot(sum.order, aes(x = Species, y = mean, fill = col_order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mycolors) +
  ylab("Relative Abundance") +
  xlab("\nHost Species") +
  ggtitle("Top Order at Ambient (T0)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#it worked! 

#can we make this more like the genus level in colour and organization? 

str(sum.order$Species)
sum.order$Species <- factor(sum.order$Species, levels = c("Montipora_capitata", "Porites_compressa", "Pavona_varians", 
                                      "Pocillopora_acuta"))
pdf(file = "../output/organised_output/6c_Order_bars_ambient_correct_color.pdf")
p.order <- ggplot(sum.order, aes(x = Species, y = mean, fill = col_order)) +
  geom_bar(stat = "identity", position = "stack") +
  #ggtitle("Communities by Phylum at Ambient Temperatures (T0)") +
  scale_fill_manual(values = c('#808080','#654522','#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032',
                               '#C2B280', '#E68FAC', '#0067A5'),  
                    labels = c("<3.4%", "Rhizobiales",
                               "Cellvibrionales", "Rickettsiales", "Alteromonadales",
                               "Flavobacteriales", "Rhodobacterales" , "Phormidesmiales", 
                               "Cytophagales", "Myxococcales","Oceanospirillales")) +
  #xlab("\nHost Species") +
  xlab("") +
  scale_x_discrete(labels= c("M. capitata", "P. compressa", "P. varians", "P. acuta")) +
  ylab("Mean Relative Abundance") +
  labs(fill = "Order") +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  #coord_flip() + #this flips the plot on its side
  theme_classic() +
  #theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.line.x = element_line(linewidth = 0.15) , axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8, face = 'bold.italic')) +
  theme(axis.line.y = element_line(linewidth = 0.15), axis.title.y = element_text(size = 12), axis.ticks.y = element_blank()) +
  theme(legend.box.margin = margin(t = 20), legend.text = element_text(size = 8))
dev.off()
# Genus level plot currently deprecated version below here
#Commented out, but could be useful in future:

#3 Can we make a class level plot??
#1 phylum level
percent.class <- bac.ambient.ra %>% 
  tax_glom(taxrank = "Class") 
perc.melt.class <- psmelt(percent.class)

sum.class <- plyr::ddply(perc.melt.class, c("Class", "Species"), summarise,
                         N = length(Abundance),
                         mean = mean(Abundance),
                         sd = sd(Abundance), 
                         se = sd/sqrt(N)
)

ggplot(sum.class, aes(x = Species, y = mean, fill = Class)) +
  geom_bar(stat = "identity") +
  ggtitle("Communities by Class at Ambient Temperatures (T0)") +
  #guides(fill=NULL) +
  xlab("/nHost Species") +
  ylab("Mean Relative Abundance") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))


x = tapply(sum.class$mean, sum.class$Class, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x) #add top 10 to code below to set NAs

##Bacteria_Unclassified is not included as top 10 below bc we can't distinguish whether they all come from the same phylum
#Added # 11 = Verrucomicrobiae, and moved unclassfied to Unassigned

sum.class$Class = factor(as.character(sum.class$Class), levels=names(x))
sum.class$col_class <- sum.class$Class
View(sum.class$col_class)
#Set everything that is super low to NA so that we can call them "other"
sum.class$col_class <- as.character(sum.class$col_class)
sum.class$col_class <- ifelse(is.na(sum.class$col_class), 
                              'Unassigned', sum.class$col_class)
sum.class$col_class <- as.factor(sum.class$col_class)

sum.class$col_class[sum.class$col_class != "Gammaproteobacteria" &
                      sum.class$col_class != "Deltaproteobacteria" & 
                      sum.class$col_class != "Bacteroidia" &
                      sum.class$col_class != "Alphaproteobacteria" &
                      sum.class$col_class != "Oxyphotobacteria" &
                      sum.class$col_class != "Campylobacteria" &
                      sum.class$col_class != "Actinobacteria" &
                      sum.class$col_class != "Spirochaetia" &
                      sum.class$col_class != "Bacilli" &
                      sum.class$col_class != "Verrucomicrobiae"] <- NA


levels(sum.class$col_class)
# add new factor
sum.class$col_class <- factor(sum.class$col_class, levels = c(levels(sum.class$col_class), "Other"))
# convert NAs to other
sum.class$col_class[is.na(sum.class$col_class)] = "Other"
sum.class$col_class <- factor(x = sum.class$col_class, levels = c("Other", "Verrucomicrobiae",
                                                                  "Bacilli", "Spirochaetia", "Actinobacteria",
                                                                  "Campylobacteria", "Oxyphotobacteria" , "Alphaproteobacteria", 
                                                                  "Bacteroidia", "Deltaproteobacteria","Gammaproteobacteria"))
#Make a colour scheme
library(RColorBrewer)
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

#plot!
ggplot(sum.class, aes(x = Species, y = mean, fill = col_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mycolors) +
  ylab("Relative Abundance") +
  xlab("\nHost Species") +
  ggtitle("Top Class at Ambient (T0)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#it worked! 

#Can we change the colours and match it to genus level plot? 
#First change order of host species
str(sum.class$Species)
sum.class$Species <- factor(sum.class$Species, levels = c("Montipora_capitata", "Porites_compressa", "Pavona_varians", 
                                                          "Pocillopora_acuta"))
pdf(file = "../output/organised_output/6d_Class_bars_ambient_correct_color.pdf")
p.class <- ggplot(sum.class, aes(x = Species, y = mean, fill = col_class)) +
  geom_bar(stat = "identity", position = "stack") +
  #ggtitle("Communities by Phylum at Ambient Temperatures (T0)") +
  scale_fill_manual(values = c('#808080','#654522','#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032',
                               '#C2B280', '#E68FAC', '#0067A5'),  
                    labels = c("<1% abundance", "Verrucomicrobiae",
                               "Bacilli", "Spirochaetia", "Actinobacteria",
                               "Campylobacteria", "Oxyphotobacteria" , "Alphaproteobacteria", 
                               "Bacteroidia", "Deltaproteobacteria","Gammaproteobacteria")) +
  #xlab("\nHost Species") +
  xlab("") +
  scale_x_discrete(labels= c("M. capitata", "P. compressa", "P. varians", "P. acuta")) +
  ylab("Mean Relative Abundance") +
  labs(fill = "Class") +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  #coord_flip() + #this flips the plot on its side
  theme_classic() +
  #theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  theme(axis.line.x = element_line(linewidth = 0.15) , axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8, face = 'bold.italic')) +
  theme(axis.line.y = element_line(linewidth = 0.15), axis.title.y = element_text(size = 12), axis.ticks.y = element_blank()) +
  theme(legend.box.margin = margin(t = 20), legend.text = element_text(size = 8))
dev.off()

##Can we put them on the same graph on top of each other? 

library(patchwork)
# Combine the plots

#Can we make the legend below??
#p.phy.amb <- p.phy.amb + theme(
#  legend.position = "bottom") +
#  guides(fill=guide_legend(ncol=2,bycol=TRUE)) + 
#  labs(fill = "")

#p.gen.amb <- p.gen.amb + theme(
#  legend.position = "bottom") +
#  #legend.box.margin = margin(b = 50)) +
#  guides(fill=guide_legend(ncol=2,bycol=TRUE)) + 
#  labs(fill = "")

#plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
gen.phy.wrap <- wrap_plots(p.phy.amb, p.gen.amb, ncol = 2, nrow = 1)
gen.class.wrap <- wrap_plots(p.class, p.gen.amb, ncol = 2, nrow = 1)
gen.order.wrap <- wrap_plots(p.order, p.gen.amb, ncol = 2, nrow = 1)
#ggsave("bac_figures/output/organised_output/5_6_bars_ambient_coord_wrap.pdf", plot = gen.phy.wrap)
ggsave("../output/organised_output/5_6_gen_class_ambient_nonflip_wrap.pdf", plot = gen.class.wrap, width = 210, height = 99, units = "mm" )
ggsave("../output/organised_output/5_6_gen_order_ambient_nonflip_wrap.pdf", plot = gen.order.wrap, width = 210, height = 99, units = "mm" )





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

###Pull out an entire Endozoicomonas spreadsheet###

rank_names(bac.exp)
endo <- subset_taxa(bac.exp, Family == "Endozoicomonadaceae")
tax_table(endo)
View(tax_table(endo))

endo.otu <- as.data.frame(otu_table(endo))
endo.data <- as.data.frame(sample_data(endo))
View(endo.data)

write.csv(endo.otu, "bac_figures/output/endozoicomonas_otu.csv")
write.csv(endo.data, "bac_figures/output/endozoicomonas_metadata.csv")

#Make and Endo tree#
endo.amb <- subset_taxa(bac.ambient, Family == "Endozoicomonadaceae")

endo.amb.mcap <- subset_samples(endo.amb, Species == "Montipora_capitata")
endo.amb.pvar <- subset_samples(endo.amb, Species == "Pavona_varians")
endo.amb.pacu <- subset_samples(endo.amb, Species == "Pocillopora_acuta")
endo.amb.pcomp <- subset_samples(endo.amb, Species == "Porites_compressa")
endo.amb
plot_tree(endo.amb.mcap, color = "Parent.ID")
p.tree <- plot_tree(endo.amb, shape = "Species", color = "Clade")
p.tree + facet_wrap(~Species)







