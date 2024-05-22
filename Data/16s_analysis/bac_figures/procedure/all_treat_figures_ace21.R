#### New R Script to streamline the 16S figure generation & analysis for ACE21 ####
### T1, TF ###

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(plyr)


##NOTE: The working colours are as follows: 
# mcap "#E69F00" ; pcomp "#56B4E9" ; pvar "#CC79A7" ; pacu"#009E73"

#E.g.,
#Colors for ambient and high:       
scale_fill_manual(values=c("deepskyblue1", "orangered"))
#Colors for the corals species:
#Mcap, Pcomp, Pvar, Pacu
scale_fill_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73"))
#Colors for these species by algal profile:
scale_fill_manual(values=c("darkorange3", "pink1","powderblue","#0072B2" ,"firebrick","darkgreen",    "palegreen2","goldenrod1"))

#Check your working directory

getwd()

#Load 16S data in & extract the samples:

load("16s_phyloseq4HE.RData")
bac.exp <- subset_samples(Bac.seq, Type == "sample")
bac.exp.data <- as(sample_data(bac.exp), "data.frame")

#Just look at richness and evenness over time, coloured by high and ambient temp?

alpha <- phyloseq::estimate_richness(bac.exp, measures = c("Observed", "Shannon"))
alpha$Evenness <- alpha$Shannon/log(alpha$Observed)
alpha$sample_name.1 <- rownames(alpha)
alpha <- merge(alpha, bac.exp.data, by = "sample_name.1")

sum.obs <- ddply(alpha, c("Species", "Time.Point", "Treatment"), summarise,
                 N = length(Observed),
                 mean = mean(Observed),
                 sd = sd(Observed), 
                 se = sd/sqrt(N)
)

pdf(file = "bac_figures/output/organised_output/7_Richness_by_Treatment.pdf")
ggplot(sum.obs, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment), size = 3) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Microbial Richness") +
  xlab("\nTime Point") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Species)
dev.off()


#Evenness
sum.even <- ddply(alpha, c("Species", "Time.Point", "Treatment"), summarise,
                 N = length(Evenness),
                 mean = mean(Evenness),
                 sd = sd(Evenness), 
                 se = sd/sqrt(N)
)


pdf(file = "bac_figures/output/organised_output/8_Evenness_by_Treatment.pdf")
ggplot(sum.even, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment), size = 3) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Microbial Evenness") +
  xlab("\nTime Point") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Species)
dev.off()


##How about extracting the distance from centroids? 
#I think we need a column that combines Treatment with Time.Point

bac.ra <- transform_sample_counts(bac.exp, function(x) x/ sum(x)) 

bac.exp.data$Treat.Time <- paste(bac.exp.data$Time.Point, bac.exp.data$Treatment, sep= "_")
bac.exp.data$Host.Treat.Time <- paste(bac.exp.data$Species, bac.exp.data$Treat.Time, sep = "_")
#Let's run a quick permanova to look at differences in bray curtis by species
bc <- phyloseq::distance(bac.ra, method = "bray")
vegan::adonis2(bc ~ Host.Treat.Time, data = bac.exp.data)
vegan::permutest(vegan::betadisper(bc, bac.exp.data$Host.Treat.Time, type = "centroid"))
disp <- vegan::betadisper(bc, bac.exp.data$Host.Treat.Time, type = "centroid")


#1) 
#Extract distances from centroid to plot in ggplot
disp.data <- as.data.frame(disp$distances)
disp.data$sample_name.1 <- rownames(disp.data)
disp.data <- merge(disp.data, bac.exp.data, by = "sample_name.1")

sum.disp <- ddply(disp.data, c("Species", "Time.Point", "Treatment"), summarise,
                  N = length(`disp$distances`),
                  mean = mean(`disp$distances`),
                  sd = sd(`disp$distances`), 
                  se = sd/sqrt(N)
)


pdf(file = "bac_figures/output/organised_output/9_Distance_to_Centroid_by_Treatment.pdf")
ggplot(sum.disp, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment), size = 3) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Distance to Centroid") +
  xlab("\nTime.Point") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Species)
dev.off()

pdf(file = "bac_figures/output/organised_output/10_Distance_to_Centroid_split_by_Time.pdf")
ggplot(disp.data, aes(x = Treatment, y = `disp$distances`)) +
  geom_boxplot(aes(color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Distance to Centroid") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_grid(Species ~ Time.Point)
dev.off()

#Okay let's make some ordinations time point
?ordinate()
ord.all.bc <- ordinate(bac.ra, "NMDS", "bray", trymax = 1000) 
#This won't converge ---- PROBLEM
vegan::stressplot(ord.all.bc)
scores(ord.all.bc, display="sites")
nmds.all.scores <- as.data.frame(ord.all.bc$sites)
nmds.all.scores$host <- bac.exp.data$Species
nmds.all.scores$symbiont <- bac.exp.data$Clade
nmds.all.scores$time.point <- bac.exp.data$Time.Point
nmds.all.scores$treatment <- bac.exp.data$Treatment
nmds.all.scores$host.treat.time <- bac.exp.data$Host.Treat.Time


hull_all_nmds <- nmds.all.scores %>%
  group_by(host.treat.time) %>%
  slice(chull(NMDS1, NMDS2))

pdf(file = "bac_figures/output/organised_output/11_NMDS_bray_all_by_time.pdf")
ggplot(nmds.all.scores, aes(x = NMDS1, y= NMDS2)) + 
  geom_point(aes(colour = nmds.all.scores$treatment), size = 3) +
  geom_polygon(data = hull_all_nmds, aes(x=NMDS1,y=NMDS2,fill=treatment),alpha=0.30) +
  scale_fill_manual(values=c( "deepskyblue1", "orangered")) +
  scale_colour_manual(values=c( "deepskyblue1", "orangered")) +
  labs(color = 'Treatment', fill = 'Treatment') +
  theme_classic() +
  #theme(legend.position = c(0.8, 0.75)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))+
  facet_grid(host ~ time.point)
dev.off()


##How about looking at top abundant taxa over time

mcap.phy <- subset_samples(bac.exp, Species == "Montipora_capitata")

pcomp.phy <- subset_samples(bac.exp, Species == "Porites_compressa")

pvar.phy <- subset_samples(bac.exp, Species == "Pavona_varians")

pacu.phy <- subset_samples(bac.exp, Species == "Pocillopora_acuta")


percent.gen <- pacu.phy %>% 
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) 
perc.melt.gen <- psmelt(percent.gen)
sum <- ddply(perc.melt.gen, c("Genus", "Species"), summarise,
                 N = length(Abundance),
                 mean = mean(Abundance),
                 sd = sd(Abundance), 
                 se = sd/sqrt(N)
)
x = tapply(sum$mean, sum$Genus, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
head(x)


#Okay what if I make 4 sets of plots

#mcap
#P3OB-42_ge 
#Endozoicomonas  
#MBIC10086 
#Prosthecochloris 
#Ruegeria
mcap.top <- subset_taxa(mcap.phy, Genus =="P3OB-42_ge" | Genus == "Endozoicomonas"| Genus == "MBIC10086" | Genus == "Prosthecochloris" | Genus =="Ruegeria")
mcap.top.ra <- transform_sample_counts(mcap.top, function(x) {x/sum(x)})
mcap.top.data <- as(sample_data(mcap.top), "data.frame")
percent.gen <- tax_glom(mcap.top.ra, taxrank = "Genus") 

perc.melt.gen <- psmelt(percent.gen)
sum <- ddply(perc.melt.gen, c("Genus", "Species", "Treatment", "Time.Point" ), summarise,
             N = length(Abundance),
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)
pdf(file = "bac_figures/output/organised_output/12a_Mcap_top5abundant_genera_bytreatment.pdf")
ggplot(data = sum, aes(x = Treatment, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_grid(Genus ~ Time.Point)
dev.off()

pdf(file = "bac_figures/output/organised_output/12b_Mcap_top5abundant_genera_bytimepoint.pdf")
ggplot(data = sum, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Genus)
dev.off()


#pcomp
#Endozoicomonas
#Candidatus_Amoebophilus
#Ruegeria
#Rhodobacteraceae_unclassified
#Francisellaceae_ge

pcomp.top <- subset_taxa(pcomp.phy, Genus == "Endozoicomonas"| Genus == "Candidatus_Amoebophilus" | Genus == "Ruegeria" | Genus == "Rhodobacteraceae_unclassified" | Genus =="Francisellaceae_ge")
pcomp.top.ra <- transform_sample_counts(pcomp.top, function(x) {x/sum(x)})
pcomp.top.data <- as(sample_data(pcomp.top), "data.frame")
percent.gen <- tax_glom(pcomp.top.ra, taxrank = "Genus") 

perc.melt.gen <- psmelt(percent.gen)
sum <- ddply(perc.melt.gen, c("Genus", "Species", "Treatment", "Time.Point" ), summarise,
             N = length(Abundance),
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)
pdf(file = "bac_figures/output/organised_output/13a_Pcomp_top5abundant_genera_bytreatment.pdf")
ggplot(data = sum, aes(x = Treatment, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_grid(Genus ~ Time.Point)
dev.off()

pdf(file = "bac_figures/output/organised_output/13b_Pcomp_top5abundant_genera_bytimepoint.pdf")
ggplot(data = sum, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Genus)
dev.off()


#pvar
#Candidatus_Amoebophilus 
#Endozoicomonas
#Cytophagales_unclassified 
#Rhodobacteraceae_unclassified 
#Fulvivirga

pvar.top <- subset_taxa(pvar.phy, Genus == "Candidatus_Amoebophilus" | Genus == "Endozoicomonas"| Genus == "Cytophagales_unclassified" | Genus == "Rhodobacteraceae_unclassified" | Genus =="Fulvivirga")
pvar.top.ra <- transform_sample_counts(pvar.top, function(x) {x/sum(x)})
pvar.top.data <- as(sample_data(pvar.top), "data.frame")
percent.gen <- tax_glom(pvar.top.ra, taxrank = "Genus") 

perc.melt.gen <- psmelt(percent.gen)
sum <- ddply(perc.melt.gen, c("Genus", "Species", "Treatment", "Time.Point" ), summarise,
             N = length(Abundance),
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)
pdf(file = "bac_figures/output/organised_output/14a_Pvar_top5abundant_genera_bytreatment.pdf")
ggplot(data = sum, aes(x = Treatment, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_grid(Genus ~ Time.Point)
dev.off()

pdf(file = "bac_figures/output/organised_output/14b_Pvar_top5abundant_genera_bytimepoint.pdf")
ggplot(data = sum, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Genus)
dev.off()


#pacu
#Candidatus_Amoebophilus  
#Endozoicomonas
#Rhodobacteraceae_unclassified 
#Algicola
#Thalassospira
pacu.top <- subset_taxa(pacu.phy, Genus == "Candidatus_Amoebophilus" | Genus == "Endozoicomonas"| Genus == "Rhodobacteraceae_unclassified" | Genus == "Algicola" | Genus =="Thalassospira")
pacu.top.ra <- transform_sample_counts(pacu.top, function(x) {x/sum(x)})
pacu.top.data <- as(sample_data(pacu.top), "data.frame")
percent.gen <- tax_glom(pacu.top.ra, taxrank = "Genus") 

perc.melt.gen <- psmelt(percent.gen)
sum <- ddply(perc.melt.gen, c("Genus", "Species", "Treatment", "Time.Point" ), summarise,
             N = length(Abundance),
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)
pdf(file = "bac_figures/output/organised_output/15a_Pacu_top5abundant_genera_bytreatment.pdf")
ggplot(data = sum, aes(x = Treatment, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_grid(Genus ~ Time.Point)
dev.off()

pdf(file = "bac_figures/output/organised_output/15b_Pacu_top5abundant_genera_bytimepoint.pdf")
ggplot(data = sum, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Genus)
dev.off()



#How about just a plot of all the Endozoicomonas

bac.ra <- transform_sample_counts(bac.exp, function(x) {x/sum(x)})
endo <- subset_taxa(bac.ra, Genus == "Endozoicomonas")
endo.data <- as(sample_data(endo), "data.frame")

percent.gen <- tax_glom(endo, taxrank = "Genus") 

perc.melt.gen <- psmelt(percent.gen)
sum <- ddply(perc.melt.gen, c("Genus", "Species", "Treatment", "Time.Point" ), summarise,
             N = length(Abundance),
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)
pdf(file = "bac_figures/output/organised_output/16a_Endo_bytreatment.pdf")
ggplot(data = sum, aes(x = Treatment, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_grid(Species ~ Time.Point)
dev.off()

pdf(file = "bac_figures/output/organised_output/16b_Endo_bytimepoint.pdf")
ggplot(data = sum, aes(x = Time.Point, y = mean)) +
  geom_point(aes(color = Treatment)) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se), width = 0.1, color = Treatment)) +
  geom_line(aes(group = Treatment, color = Treatment)) +
  scale_color_manual(values=c( "deepskyblue1", "orangered")) +
  theme_classic() +
  ylab("Mean Relative Abundance") +
  xlab("\nTreatment") +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap(~Species)
dev.off()
