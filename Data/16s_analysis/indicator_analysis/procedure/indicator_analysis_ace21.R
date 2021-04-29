###Multi-level pattern association analysis for Matsuda et al. 2021###
##Indicator species analysis##


library(indicspecies)
library(phyloseq)
library(plyr)

setwd("~/Documents/OSUDocs/Projects/ACE21/Hawaiian-coral-thermal-tolerance/Data/16s_analysis/")
load("16s_phyloseq4HE.RData")

##Split data by treatment

bac.exp <- subset_samples(Bac.seq, Type == "sample")

#The next steps are kinda useless unless we want to test t0 as ambient
#Concatenate & make a new data column to combine all data at T0 as "Ambient" 
sample_data(bac.exp)['Treat.Time'] <- paste(sample_data(bac.exp)$Treatment, sample_data(bac.exp)$Time.Point)
View(as(sample_data(bac.exp), "data.frame"))
#Revalue
sample_data(bac.exp)$Treat.Time <- as.factor(sample_data(bac.exp)$Treat.Time)
levels(sample_data(bac.exp)$Treat.Time)
sample_data(bac.exp)$Treat.Time <- revalue(sample_data(bac.exp)$Treat.Time, c("Ambient T0" = "Ambient", "Ambient T1" = "Ambient", "Ambient TF" = "Ambient", "High T0" = "Ambient", 
                                                                              "High T1" = "High", "High TF" = "High"))
levels(sample_data(bac.exp)$Treat.Time)

##Let's try it on just the T0 dataset between species
t0 <- subset_samples(bac.exp, Treatment == "T0")
t0 <- prune_taxa(taxa_sums(t0) > 0, t0)
t0.otu <- as(otu_table(t0), "matrix")
if(taxa_are_rows(t0)) {t0.otu <- t(t0.otu)}
t0.otu <- as.data.frame(t0.otu)
t0.data <- as(sample_data(t0), "data.frame")

indval.t0 <- multipatt(t0.otu, t0.data$Species, control = how(nperm = 999))
#view results
summary(indval.t0)

###Here I copied the summary to a .csv file & re-uploaded here as indval.t0.sum
indval.t0.sum <- read.csv("indicator_analysis/output/species_indicval.csv")

##Add taxonomy
tax<-as(tax_table(t0),"matrix")
head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
head(tax)

#Merge taxonomy with your indicator values
indval.t0.sum.tax <- merge(indval.t0.sum,tax, by = "taxa_id")

head(indval.t0.sum.tax)

write.csv(indval.t0.sum.tax, "indicator_analysis/output/species_indicval_with_tax.csv")
