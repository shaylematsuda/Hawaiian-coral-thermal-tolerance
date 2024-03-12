##Creating the Figure of Top Abundant + Top Core + Top Differential by Host species
##For Matsuda et al. 2021

#load libraries
library(phyloseq)
library(ggplot2)
library(tidyr)

#Due diligence: Make sure you're in the procedure directory
getwd()

#bring in data from ANCOM
mcap.ancom.melt <- read.csv("../../ancom/output/by_treatment/mcap_ancom.csv")
pcomp.ancom.melt <- read.csv("../../ancom/output/by_treatment/pcomp_ancom.csv")
pvar.ancom.melt <- read.csv("") #So there are no significantly differentially abundant taxa for p var by treatment
pacu.ancom.melt <- read.csv("../../ancom/output/by_treatment/pacu_ancom.csv")

#bring in data from CORE
mcap.core.melt <- read.csv("../../core_microbiome/output/mcap_core.csv")
pcomp.core.melt <- read.csv("../../core_microbiome/output/pcomp_core.csv")
pvar.core.melt <- read.csv("../../core_microbiome/output/pvar_core.csv")
pacu.core.melt <- read.csv("../../core_microbiome/output/pacu_core.csv")

#bring in the data from top abundant taxa
mcap.abund.melt <- read.csv("../../top_abundance/output/mcap_abund.csv")
pcomp.abund.melt <-read.csv("../../top_abundance/output/pcomp_abund.csv")
pvar.abund.melt <- read.csv("../../top_abundance/output/pvar_abund.csv")
pacu.abund.melt <- read.csv("../../top_abundance/output/pacu_abund.csv")


##1 Montipora Capitata
#First this is a work around to get a list of the samples in the right order for plotting 
#in order of: time, treatment, clade
mcap.toget <- mcap.core.melt[order(mcap.core.melt$Time.Point, mcap.core.melt$Treatment, mcap.core.melt$Clade),]
mcap.toget <- filter(mcap.toget, OTU == "Otu00001")
mcap.list <- as.list(mcap.toget$Sample)

#merge all data
mcap.top <- rbind(mcap.ancom.melt, mcap.core.melt, mcap.abund.melt)
mcap.top <- mcap.top %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(mcap.top, "bac_figures/output/mcap_top_taxa.csv")


#OG bubble plots
p.mcap <- ggplot(mcap.top, aes(x = Sample, y = otu_genus, size = Abundance, color = Time.Point, shape = Treatment)) +
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


##New plot trial based on treatment 
#and what if we remove the differential bc its kinda hard to read??
library(plyr)
mcap.top <- rbind(mcap.core.melt, mcap.abund.melt)
mcap.top <- rbind(mcap.ancom.melt, mcap.core.melt, mcap.abund.melt)
mcap.top <- mcap.top %>% unite(OTU, Genus, col='otu_genus',sep='-')

sum.mcap <- ddply(mcap.top, c("Treatment", "top.group", "Time.Point", "otu_genus"), summarise,
                  N = length(Abundance), 
                  mean = mean(Abundance),
                  sd = sd(Abundance), 
                  se = sd/sqrt(N)
)

ggplot(sum.mcap, aes(x = Treatment, y = mean, group = otu_genus, color = otu_genus)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_errorbar(aes(ymin = (mean-se), ymax = (mean+se)), alpha = .5, width = 0)+
  geom_line(alpha = 0.5)+
  ylab("Mean Relative Abundance \n") +
  xlab("\n Treatment") +
  labs(color='OTU Genus') +
  #scale_color_discrete(breaks = c("Otu00001-Endozoicomonas", "Otu00002-Endozoicomonas", "Otu00031-Endozoicomonas",
                                 #"Otu00078-Endozoicomonas", "Otu00004-P3OB-42_ge","Otu00013-P3OB-42_ge","Otu00015-MBIC10086",
                                 #"Otu00025-Ruegeria", "Otu00028-Candidatus_Amoebophilus", "Otu00058-Prosthecochloris", "Otu00092-uncultured",
                                 #"Otu00139-uncultured")) +
  #scale_color_manual(labels = c("Otu00001-Endozoicomonas", "Otu00002-Endozoicomonas", "Otu00031-Endozoicomonas",
                                #"Otu00078-Endozoicomonas", "Otu00004-P3OB-42_ge","Otu00013-P3OB-42_ge","Otu00015-MBIC10086",
                                #"Otu00025-Ruegeria", "Otu00028-Candidatus_Amoebophilus", "Otu00058-Prosthecochloris", "Otu00092-uncultured",
                                #"Otu00139-uncultured"),
                     #values = c("forestgreen", "green2", "green", "limegreen", "purple", "violet", "orange", "brown", "yellow", "grey", "blue3", "blue")) +
  theme_bw() +
  #facet_wrap(~Time.Point)
  facet_grid(top.group ~ Time.Point)

##what if we join EVERYTHING together and then facet grid by species? 
mcap.top <- rbind(mcap.ancom.melt, mcap.core.melt, mcap.abund.melt)
pcomp.top <- rbind(pcomp.ancom.melt, pcomp.core.melt, pcomp.abund.melt)
pvar.top <- rbind(pvar.core.melt, pvar.abund.melt) #no ancom results
pacu.top <- rbind(pacu.ancom.melt, pacu.core.melt, pacu.abund.melt) #check below code for overwrite of the error

all.top <- rbind(mcap.top, pcomp.top, pvar.top)
all.top <- all.top[,-1]
all.top <- rbind(all.top, pacu.top)

#add an otu_genus column
all.top <- all.top %>% unite(OTU, Genus, col='otu_genus',sep='-')

sum <- ddply(all.top, c("Species", "top.group", "Treatment", "Time.Point", "otu_genus"), summarise,
                  N = length(Abundance), 
                  mean = mean(Abundance),
                  sd = sd(Abundance), 
                  se = sd/sqrt(N)
)

ggplot(sum, aes(x = Treatment, y = mean, group = otu_genus, color = otu_genus)) +
  geom_point(alpha = 0.5, size = 3, aes(shape = top.group)) +
  geom_errorbar(aes(ymin = (mean-se), ymax = (mean+se)), alpha = .5, width = 0)+
  geom_line(alpha = 0.5)+
  ylab("Mean Relative Abundance \n") +
  xlab("\n Treatment") +
  labs(color='OTU Genus') +
  theme_bw() +
  facet_wrap(~Time.Point) +
  facet_grid(Species ~ Time.Point)
ggsave(filename = "../output/all_treatment_by_abund.pdf", plot = last_plot())
dev.off()



##Can we split between abundant & core? 
#OR do each of the three for each time point? 




##Thoughts - can we group the OTUs of the same genus together so that colours make more sense?
#Hang tight on this thought

#let's make a venn diagram
##Venn diagram at T0 only - how many of these taxa overlap? 
##Just the ones we pulled - interesting at OTU level AND genus level/family level
#install.packages("VennDiagram") #if needed install
library(VennDiagram)

#Need to make 4 lists to run a venn diagram

mcap.list <- mcap.top$otu_genus
pcomp.list <- all.top %>% filter(Species == "Porites_compressa") 
pcomp.list <- pcomp.list$otu_genus
pvar.list <- all.top %>% filter(Species == "Pavona_varians")
pvar.list <- pvar.list$otu_genus
pacu.list <- all.top %>% filter(Species == "Pocillopora_acuta")
pacu.list <- pacu.list$otu_genus
candidates <- list("Montipora capitata" = mcap.list, "Porites compressa" = pcomp.list,
                   "Pavona varians" = pvar.list, "Pocillopora acuta" = pacu.list)
grid.draw(venn.diagram(x = candidates, filename = NULL, fill = c("#E69F00", "#56B4E9", "#CC79A7","#009E73"), alpha = 0.5))

#Okay this works! Can we do 3 separate ones for TO only, that we can separate by top group?

all.top.T0 <- all.top %>% filter(Time.Point == "T0")

#1. Top 10 most abundant taxa
all.top.abund <- all.top.T0 %>% filter(top.group == "abund")

mcap.list.abund <- all.top.abund %>% filter(Species == "Montipora_capitata")
mcap.list.abund <- unique(mcap.list.abund$otu_genus)
pcomp.list.abund <- all.top.abund %>% filter(Species == "Porites_compressa")
pcomp.list.abund <- unique(pcomp.list.abund$otu_genus)
pvar.list.abund <- all.top.abund %>% filter(Species == "Pavona_varians")
pvar.list.abund <- unique(pvar.list.abund$otu_genus)
pacu.list.abund <- all.top.abund %>% filter(Species == "Pocillopora_acuta")
pacu.list.abund <- unique(pacu.list.abund$otu_genus)
candidates.abund <- list("Montipora capitata" = mcap.list.abund, "Porites compressa" = pcomp.list.abund,
                         "Pavona varians" = pvar.list.abund, "Pocillopora acuta" = pacu.list.abund)
overlap <- calculate.overlap(candidates.abund)
names(overlap) <- c("a1234", "a123", "a124", "a134", "a234", "a12", "a13", "a14", "a23", "a24", "a34", "a1", "a2", "a3", "a4")
View(overlap)
#Still callenging to know which quandrant is which. 1, 2, 3, 4 refers to each group but unclear which
#group they are... 
#a1234 is the one that overlaps in all of them. a123 is all that overlap in 3 of the groups (minus the one overlapping in all)

dev.off()
venn <- venn.diagram(x = candidates.abund, filename = NULL, fill = c("#E69F00", "#56B4E9", "#CC79A7","#009E73"), alpha = 0.5)
grid.draw(venn)
dev.off()
#To save
pdf(file = "../output/venn_topabundant_T0.pdf")
grid.draw(venn)
dev.off()

#2. Core
all.core <- all.top.T0 %>% filter(top.group == "core")

mcap.list.core <- all.core %>% filter(Species == "Montipora_capitata")
mcap.list.core <- unique(mcap.list.core$otu_genus)
pcomp.list.core <- all.core %>% filter(Species == "Porites_compressa")
pcomp.list.core <- unique(pcomp.list.core$otu_genus)
pvar.list.core <- all.core %>% filter(Species == "Pavona_varians")
pvar.list.core <- unique(pvar.list.core$otu_genus)
pacu.list.core <- all.core %>% filter(Species == "Pocillopora_acuta")
pacu.list.core <- unique(pacu.list.core$otu_genus)

candidates.core <- list("Montipora capitata" = mcap.list.core, "Porites compressa" = pcomp.list.core,
                         "Pavona varians" = pvar.list.core, "Pocillopora acuta" = pacu.list.core)
overlap <- calculate.overlap(candidates.core)
names(overlap) <- c("a1234", "a123", "a124", "a134", "a234", "a12", "a13", "a14", "a23", "a24", "a34", "a1", "a2", "a3", "a4")
View(overlap)
##the numbers 1,2,3,4 refer to each of the four lists in order
#meaning that 1 = M. cap, 2 = P. comp, 3 = P. var and 4 = P. acu
#e.g., "a1234" = all four together, "a134" = M cap, P. var and P. acu

venn <- venn.diagram(x = candidates.core, filename = NULL, fill = c("#E69F00", "#56B4E9", "#CC79A7","#009E73"), alpha = 0.5)
grid.draw(venn)
dev.off()
#TO save
pdf(file = "../output/venn_core_T0.pdf")
grid.draw(venn)
dev.off()

#3. Differential @ T1??

##OKAY something weird here - the way I brought in the data means that they don't differentiate by time point
##Need to go back to original abundant vs. core vs. differential to fix this. 
#UGH

all.top.T1 <- all.top %>% filter(Time.Point == "T1")
all.dif <- all.top.T1 %>% filter(top.group == "differential")

mcap.list.dif <- all.dif %>% filter(Species == "Montipora_capitata")
mcap.list.dif <- unique(mcap.list.dif$otu_genus)
pcomp.list.dif <- all.dif %>% filter(Species == "Porites_compressa")
pcomp.list.dif <- unique(pcomp.list.dif$otu_genus)
pvar.list.dif <- all.dif %>% filter(Species == "Pavona_varians")
pvar.list.dif <- unique(pvar.list.dif$otu_genus)
pacu.list.dif <- all.dif %>% filter(Species == "Pocillopora_acuta")
pacu.list.dif <- unique(pacu.list.dif$otu_genus)

candidates.dif <- list("Montipora capitata" = mcap.list.dif, "Porites compressa" = pcomp.list.dif,
                        "Pavona varians" = pvar.list.dif, "Pocillopora acuta" = pacu.list.dif)
overlap <- calculate.overlap(candidates.dif)
names(overlap) <- c("a1234", "a123", "a124", "a134", "a234", "a12", "a13", "a14", "a23", "a24", "a34", "a1", "a2", "a3", "a4")
View(overlap)
#Still callenging to know which quandrant is which. 1, 2, 3, 4 refers to each group but unclear which
#group they are... 
#a1234 is the one that overlaps in all of them. a123 is all that overlap in 3 of the groups (minus the one overlapping in all)

venn <- venn.diagram(x = candidates.dif, filename = NULL, fill = c("#E69F00", "#56B4E9", "#CC79A7","#009E73"), alpha = 0.5)
grid.draw(venn)
dev.off()
#TO save
pdf(file = "../output/venn_core_T0.pdf")
grid.draw(venn)
dev.off()



#######
##Colors to use via Shayle::

#Colors for ambient and high:       
scale_fill_manual(values=c("deepskyblue1", "orangered"))
#Colors for the corals species:
#Mcap, Pcomp, Pvar, Pacu
scale_fill_manual(values=c("#E69F00", "#56B4E9", "#CC79A7","#009E73"))
#Colors for these species by algal profile:
scale_fill_manual(values=c("darkorange3", "pink1","powderblue","#0072B2" ,"firebrick","darkgreen",    "palegreen2","goldenrod1"))

##changing the species names to have no _
mutate(Species = str_replace(Species, "_", " "))
#changing the species so they always come up in order 
meta$Species <- factor(meta$Species, levels = c("Montipora capitata", "Porites compressa", "Pavona varians","Pocillopora acuta"))
#######

##2 Porites compressa
#First this is a work around to get a list of the samples in the right order for plotting 
#in order of: time, treatment, clade
pcomp.toget <- pcomp.ancom.melt[order(pcomp.ancom.melt$Time.Point, pcomp.ancom.melt$Treatment, pcomp.ancom.melt$Clade),]
pcomp.toget <- filter(pcomp.toget, OTU == "Otu00002")
pcomp.list <- as.list(pcomp.toget$Sample)

#merge all data
pcomp.top <- rbind(pcomp.ancom.melt, pcomp.core.melt, pcomp.abund.melt)
pcomp.top <- pcomp.top %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(pcomp.top, "bac_figures/output/pcomp_top_taxa.csv")

p.pcomp <- ggplot(pcomp.top, aes(x = Sample, y = otu_genus, size = Abundance, color = Time.Point, shape = Treatment)) +
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
pvar.top <- pvar.top %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(pvar.top, "bac_figures/output/pvar_top_taxa.csv")

p.pvar<- ggplot(pvar.top, aes(x = Sample, y = otu_genus, size = Abundance, color = Time.Point, shape = Treatment)) +
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
pacu.toget <- pacu.abund.melt[order(pacu.abund.melt$Time.Point, pacu.abund.melt$Treatment, pacu.abund.melt$Clade),]
pacu.toget <- filter(pacu.toget, OTU == "Otu00001")
pacu.list <- as.list(pacu.toget$Sample)

#merge all data
pacu.ancom.melt$X <- rownames(pacu.ancom.melt) #apparently ancom doesn't have an "X" variable (but it's a dummy variable so it's fine to add or remove)
pacu.abund.melt$X <- rownames(pacu.abund.melt)
pacu.ancom.melt$X <- NULL
pacu.ancom.melt$X.1 <- NULL
pacu.abund.melt$X <- NULL
pacu.core.melt$X <- NULL
pacu.top <- rbind(pacu.ancom.melt, pacu.core.melt, pacu.abund.melt)
##NA generated..but I can't find it in the dataframe???

#Add a column to label by OTU and by genus so we can have taxonomic annotation
pacu.top <- pacu.top %>% unite(OTU, Genus, col='otu_genus',sep='-')
write.csv(pacu.top, "bac_figures/output/pacu_top_taxa.csv")

p.pacu <- ggplot(pacu.top, aes(x = Sample, y = otu_genus, size = Abundance, color = Time.Point, shape = Treatment)) +
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

