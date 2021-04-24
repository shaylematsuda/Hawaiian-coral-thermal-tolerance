###Code for ANCOM: Analysis of Composition of Microbiomes for Matsuda et al. 2021###
##Purpose: to look for taxa driving differences in 16S communities between high temp & ambient treatments

#Load libraries
library(phyloseq)
library(plyr)

#Set working directory & bring in RData
#source the ANCOM functions R script
setwd("~/Documents/OSUDocs/Projects/ACE21/Hawaiian-coral-thermal-tolerance/Data/16s_analysis/")
load("16s_phyloseq4HE.RData")
source("ancom/procedure/ANCOM_functions.R")

##Set a colour scheme for the whole script
#I kinda like these kelly colours from Becca Maher/Grace Klinges
kelly_colors = c('#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032', 
                 '#C2B280',  '#222222','#848482',  '#E68FAC', '#0067A5', 
                 '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', 
                 '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26', 
                 '#222222', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', 
                 '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5')

##Pull out montipora Capitata and porites compressa (only two host species that showed sig differences between treatment)
#Montipora capitata
mcap <- subset_samples(Bac.seq, Species == "Montipora_capitata")
mcap.t1 <- subset_samples(mcap, Time.Point == "T1")
mcap.tf <- subset_samples(mcap, Time.Point == "TF")
#Remove any taxa that have 0s across all samples
mcap <- prune_taxa(taxa_sums(mcap) > 0, mcap)
mcap.t1 <- prune_taxa(taxa_sums(mcap.t1) > 0, mcap.t1)
mcap.tf <- prune_taxa(taxa_sums(mcap.tf) > 0, mcap.tf)

#Porites compressa
pcomp <- subset_samples(Bac.seq, Species == "Porites_compressa")
pcomp.t1 <- subset_samples(pcomp, Time.Point == "T1")
pcomp.tf <- subset_samples(pcomp, Time.Point == "TF")
#Remove any taxa that have 0s across all samples
pcomp <- prune_taxa(taxa_sums(pcomp) > 0, pcomp)
pcomp.t1 <- prune_taxa(taxa_sums(pcomp.t1) > 0, pcomp.t1)
pcomp.tf <- prune_taxa(taxa_sums(pcomp.tf) > 0, pcomp.tf)

#Export OTU tables & transpose so taxa are rows
OTU.mcap <- t(otu_table(mcap))
OTU.mcap.t1 <- t(otu_table(mcap.t1))
OTU.mcap.tf <- t(otu_table(mcap.tf))

OTU.pcomp <- t(otu_table(pcomp))
OTU.pcomp.t1 <- t(otu_table(pcomp.t1))
OTU.pcomt.tf <- t(otu_table(pcomp.tf))


###ANCOM pipeline###
##Step 1: set all your parameters
#I just run through these steps for each of the time points (run two separate ANCOMs)
#then I combine downstream to make a faceted graph
#M cap T1
feature_table <- OTU.mcap.t1
meta_data <- sample_data(mcap.t1)
#M cap TF
feature_table <- OTU.mcap.tf
meta_data <- sample_data(mcap.tf)


levels(meta_data$Treatment) #Ambient first, High Second = High will be positive 
meta_data$Sample.ID <- as.character(meta_data$sample_name.1) #make sure your sample id is chr
sample_var <- "Sample.ID"
group_var <- NULL
out_cut <- 0.05  #numerical fraction, below 5% = outlier zeros, above 95% = outlier values
zero_cut <- 0.90 #numerical fraction, taxa with proportion of zeros above 90% are removed
lib_cut <- 1000 #number, removes any samples with less than lib_cut reads
neg_lb <- TRUE

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var,
                                     group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM

main_var <- "Treatment"
p_adj_method <- "BH" #default: Benjamini-Hochberg procedure
alpha <- 0.05 #level of significance
adj_formula <- NULL
rand_formula <- NULL #add in random effect 
#The random effect is a problem when split up by time point. Is there a way to split after the ANCOM? 

res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)

#Extract your data frames 
resdf.mcap.t1 <- as.data.frame(res$out)
resdf.mcap.tf <- as.data.frame(res$out)

#rbind together
resdf.mcap <- rbind(resdf.mcap.t1, resdf.mcap.tf)

# compiling results from each contrast for the figure
figdf.mcap.t1 <- as.data.frame(res$fig$data)
figdf.mcap.t1$contrast <- "T1"
rownames(figdf.mcap.t1) <- NULL

figdf.mcap.tf <- as.data.frame(res$fig$data)
figdf.mcap.tf$contrast <- "TF"
rownames(figdf.mcap.tf) <- NULL

#rbind your two datasets together
figdf.mcap <- rbind(figdf.mcap.t1, figdf.mcap.tf)

# add taxonomy
tax<-as(tax_table(mcap),"matrix")
head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
head(tax)

#Merge taxonomy with your figure data
figdf.mcap <- merge(figdf.mcap,tax, by = "taxa_id")
resdf.mcap <- merge(resdf.mcap, tax, by = "taxa_id")

head(resdf.mcap)
head(figdf.mcap)
write.csv(resdf.mcap, "ancom/output/ancom_mcap_res.csv")

# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(figdf.mcap), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(figdf.mcap$x), y = cut_off["detected_0.6"], label = "W[0.6]")

##Specialised Plot
# order genus
x = tapply(figdf.mcap$y, figdf.mcap$Genus, function(x) max(x))
x = sort(x, TRUE)
figdf.mcap$Genus = factor(as.character(figdf.mcap$Genus), levels=names(x))
figdf.mcap$col_genus <- figdf.mcap$Genus

#Set everything that is super low to NA so that we can call them "other"
figdf.mcap$col_genus[figdf.mcap$col_genus != "Alteromonas" & 
                        figdf.mcap$col_genus != "Pseudoalteromonas" &
                        figdf.mcap$col_genus != "Ruegeria" &
                        figdf.mcap$col_genus != "Reichenbachiella" &
                        figdf.mcap$col_genus != "Alteromonadales_unclassified" &
                        figdf.mcap$col_genus != "uncultured" &
                        figdf.mcap$col_genus != "Prosthecochloris" &
                        figdf.mcap$col_genus != "JGI_0000069-P22_ge"] <- NA


levels(figdf.mcap$col_genus)
# add new factor
figdf.mcap$col_genus <- factor(figdf.mcap$col_genus, levels = c(levels(figdf.mcap$col_genus), "Other"))
# convert NAs to other
figdf.mcap$col_genus[is.na(figdf.mcap$col_genus)] = "Other"

mcap.p <- ggplot(figdf.mcap, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 3) +
  facet_grid(~contrast) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  #geom_hline(yintercept = c(18, linetype = "dashed") + 
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") + #at 330.6
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 1, color = "orange", parse = TRUE) +
  ggtitle("Montipora capitata")
ggsave("ancom/output/mcap_ancom.pdf")





##All the same steps for Porites compressa! 

##Step 1: set all your parameters
#I just run through these steps for each of the time points
#then I combine downstream to make a faceted graph
#P comp T1
feature_table <- OTU.pcomp.t1
meta_data <- sample_data(pcomp.t1)
#P comp TF
feature_table <- OTU.pcomt.tf
meta_data <- sample_data(pcomp.tf)


levels(meta_data$Treatment) #Ambient first, High Second = High will be positive 
meta_data$Sample.ID <- as.character(meta_data$sample_name.1) #make sure your sample id is chr
sample_var <- "Sample.ID"
group_var <- NULL
out_cut <- 0.05  #numerical fraction, below 5% = outlier zeros, above 95% = outlier values
zero_cut <- 0.90 #numerical fraction, taxa with proportion of zeros above 90% are removed
lib_cut <- 1000 #number, removes any samples with less than lib_cut reads
neg_lb <- TRUE

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var,
                                     group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM

main_var <- "Treatment"
p_adj_method <- "BH" #default: Benjamini-Hochberg procedure
alpha <- 0.05 #level of significance
adj_formula <- NULL
rand_formula <- NULL #add in random effect 
#The random effect is a problem when split up by time point. Is there a way to split after the ANCOM? 

res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)


resdf.pcomp.t1 <- as.data.frame(res$out)
resdf.pcomp.tf <- as.data.frame(res$out)

#rbind together
resdf.pcomp <- rbind(resdf.pcomp.t1, resdf.pcomp.tf)

# compiling results from each contrast for the figure
figdf.pcomp.t1 <- as.data.frame(res$fig$data)
figdf.pcomp.t1$contrast <- "T1"
rownames(figdf.pcomp.t1) <- NULL

figdf.pcomp.tf <- as.data.frame(res$fig$data)
figdf.pcomp.tf$contrast <- "TF"
rownames(figdf.pcomp.tf) <- NULL

#rbind your two datasets together
figdf.pcomp<- rbind(figdf.pcomp.t1, figdf.pcomp.tf)

# add taxonomy
tax<-as(tax_table(pcomp),"matrix")
head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
head(tax)

#Merge taxonomy with your figure data
figdf.pcomp <- merge(figdf.pcomp,tax, by = "taxa_id")
resdf.pcomp <- merge(resdf.pcomp, tax, by = "taxa_id")

head(resdf.pcomp)
head(figdf.pcomp)
write.csv(resdf.pcomp, "ancom/output/ancom_pcomp_res.csv")

# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(figdf.pcomp), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(figdf.pcomp$x), y = cut_off["detected_0.6"], label = "W[0.6]")

##Specialised Plot
# order genus
x = tapply(figdf.pcomp$y, figdf.pcomp$Genus, function(x) max(x))
x = sort(x, TRUE)
figdf.pcomp$Genus = factor(as.character(figdf.pcomp$Genus), levels=names(x))
figdf.pcomp$col_genus <- figdf.pcomp$Genus

figdf.pcomp$col_genus[figdf.pcomp$col_genus != "Francisellaceae_ge" & 
                       figdf.pcomp$col_genus != "Francisellaceae_unclassified" &
                       figdf.pcomp$col_genus != "Endozoicomonas" &
                       figdf.pcomp$col_genus != "Proteobacteria_unclassified" &
                       figdf.pcomp$col_genus != "Rhodobacteraceae_unclassified" &
                       figdf.pcomp$col_genus != "Bythopirellula" ] <- NA


levels(figdf.pcomp$col_genus)
# add new factor
figdf.pcomp$col_genus <- factor(figdf.pcomp$col_genus, levels = c(levels(figdf.pcomp$col_genus), "Other"))
# convert NAs to other
figdf.pcomp$col_genus[is.na(figdf.pcomp$col_genus)] = "Other"

pcomp.p <- ggplot(figdf.pcomp, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 3) +
  facet_grid(~contrast) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 1, color = "orange", parse = TRUE) +
  ggtitle("Porites compressa")
ggsave("ancom/output/pcomp_ancom.pdf")


#############
##What about contrasts for M cap - C vs CD

#subset to ambient only
mcap.amb <- subset_samples(mcap, Treatment == "Ambient")
#Remove any taxa that have 0s across all samples
mcap.amb <- prune_taxa(taxa_sums(mcap.amb) > 0, mcap.amb)
OTU.mcap.amb <- t(otu_table(mcap.amb))

##Step 1: set all your parameters
#I just run through these steps for each of the time points (run two separate ANCOMs)
#then I combine downstream to make a faceted graph

#M cap Ambient only 
feature_table <- OTU.mcap.amb
meta_data <- sample_data(mcap.amb)

levels(meta_data$Clade) #C first, CD Second = CD will be positive 
meta_data$Sample.ID <- as.character(meta_data$sample_name.1) #make sure your sample id is chr
sample_var <- "Sample.ID"
group_var <- NULL
out_cut <- 0.05  #numerical fraction, below 5% = outlier zeros, above 95% = outlier values
zero_cut <- 0.90 #numerical fraction, taxa with proportion of zeros above 90% are removed
lib_cut <- 1000 #number, removes any samples with less than lib_cut reads
neg_lb <- TRUE

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var,
                                     group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM

main_var <- "Clade"
p_adj_method <- "BH" #default: Benjamini-Hochberg procedure
alpha <- 0.05 #level of significance
adj_formula <- "Time.Point"
rand_formula <- NULL #add in random effect 
#The random effect is a problem...why?

res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)

#Extract your data frame 
resdf.mcap.amb <- as.data.frame(res$out)
rownames(resdf.mcap.amb) <- NULL

# compiling results from each contrast for the figure
figdf.mcap.amb <- as.data.frame(res$fig$data)
rownames(figdf.mcap.amb) <- NULL

# add taxonomy
tax<-as(tax_table(mcap.amb),"matrix")
head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
head(tax)

#Merge taxonomy with your figure data
figdf.mcap.amb <- merge(figdf.mcap.amb,tax, by = "taxa_id")
resdf.mcap.amb <- merge(resdf.mcap.amb, tax, by = "taxa_id")

head(resdf.mcap)
head(figdf.mcap)
write.csv(resdf.mcap.amb, "ancom/output/ancom_mcap_ambclad_res.csv")

# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(figdf.mcap.amb$x), y = cut_off["detected_0.6"], label = "W[0.6]")

##Specialised Plot
# order genus
x = tapply(figdf.mcap.amb$y, figdf.mcap.amb$Genus, function(x) max(x))
x = sort(x, TRUE)
figdf.mcap.amb$Genus = factor(as.character(figdf.mcap.amb$Genus), levels=names(x))
figdf.mcap.amb$col_genus <- figdf.mcap.amb$Genus

#Set everything that is super low to NA so that we can call them "other"
figdf.mcap.amb$col_genus[figdf.mcap.amb$col_genus != "Endozoicomonas" & 
                       figdf.mcap.amb$col_genus != "Pseudoalteromonas" &
                       figdf.mcap.amb$col_genus != "Oxyphotobacteria_unclassified" ] <- NA


levels(figdf.mcap.amb$col_genus)
# add new factor
figdf.mcap.amb$col_genus <- factor(figdf.mcap.amb$col_genus, levels = c(levels(figdf.mcap.amb$col_genus), "Other"))
# convert NAs to other
figdf.mcap.amb$col_genus[is.na(figdf.mcap.amb$col_genus)] = "Other"

mcap.amb.p <- ggplot(figdf.mcap.amb, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 3) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 1, color = "orange", parse = TRUE) +
  ggtitle("Montipora capitata by clade")
ggsave("ancom/output/mcap_ambclade_ancom.pdf")


##What about a plot of endo? Just because I am curious about the mcap ambient results
mcap.amb.ra <- transform_sample_counts(mcap.amb, function(x) x/ sum(x))
endo <- subset_taxa(mcap.amb.ra, Genus == "Endozoicomonas")
merged.endo <- tax_glom(endo, "Genus", NArm = FALSE)
endo.melt <- psmelt(merged.endo)
sum.endo <- ddply(endo.melt, c("Time.Point",  "Clade"), summarise,
                     N = length(Abundance), 
                     mean = mean(Abundance),
                     sd = sd(Abundance), 
                     se = sd/sqrt(N)
)
ggplot(sum.endo, aes(x = Clade, y = mean, fill = Time.Point)) +
  geom_point(aes(color = Time.Point)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = Time.Point), width = 0.1) +
  facet_wrap(~Time.Point)
ggsave("ancom/output/mcap_ambclade_endo.pdf")

##Pavona -  C1 vs. C27 
##What about pavona vs acuta 
