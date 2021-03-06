###Code for ANCOM: Analysis of Composition of Microbiomes for Matsuda et al. 2021###
##Purpose: to look for taxa driving differences in 16S communities between high temp & ambient treatments

#Load libraries
library(phyloseq)
library(plyr)

#Set working directory & bring in RData
#source the ANCOM functions R script
setwd("~/Documents/OSUDocs/Projects/ACE21/Hawaiian-coral-thermal-tolerance/Data/16s_analysis/")
load("16s_phyloseq4HE.RData")
source("ancom/procedure/ANCOM_v2.1.R")

##Set a colour scheme for the whole script
#I kinda like these kelly colours from Becca Maher/Grace Klinges
kelly_colors = c('#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032', 
                 '#C2B280',  '#222222','#848482',  '#E68FAC', '#0067A5', 
                 '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', 
                 '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26')


##Pull out data from each species
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

#Pavona varians
pvar <- subset_samples(Bac.seq, Species == "Pavona_varians")
pvar.t1 <- subset_samples(pvar, Time.Point == "T1")
pvar.tf <- subset_samples(pvar, Time.Point == "TF")
pvar <- prune_taxa(taxa_sums(pvar) > 0, pvar)
pvar.t1 <- prune_taxa(taxa_sums(pvar.t1) > 0, pvar.t1)
pvar.tf <- prune_taxa(taxa_sums(pvar.tf) > 0, pvar.tf)

#Pocillopora acuta
#There is one poc sample that is labelled as "F1" which should be "T1" - oops!
pacu <- subset_samples(Bac.seq, Species == "Pocillopora_acuta")
levels(sample_data(pacu)$Time.Point)[levels(sample_data(pacu)$Time.Point)=="F1"] <- "T1"
pacu.t1 <- subset_samples(pacu, Time.Point == "T1")
pacu.tf <- subset_samples(pacu, Time.Point == "TF")
pacu <- prune_taxa(taxa_sums(pacu) > 0, pacu)
pacu.t1 <- prune_taxa(taxa_sums(pacu.t1) > 0, pacu.t1)
pacu.tf <- prune_taxa(taxa_sums(pacu.tf) > 0, pacu.tf)

#Export OTU tables & transpose so taxa are rows
#Monitpora capitata
OTU.mcap <- t(otu_table(mcap))
OTU.mcap.t1 <- t(otu_table(mcap.t1))
OTU.mcap.tf <- t(otu_table(mcap.tf))
#Porites compressa
OTU.pcomp <- t(otu_table(pcomp))
OTU.pcomp.t1 <- t(otu_table(pcomp.t1))
OTU.pcomp.tf <- t(otu_table(pcomp.tf))
#Pavona varians
OTU.pvar <- t(otu_table(pvar))
OTU.pvar.t1 <- t(otu_table(pvar.t1))
OTU.pvar.tf <- t(otu_table(pvar.tf))
#Pocillopora acuta
OTU.pacu <- t(otu_table(pacu))
OTU.pacu.t1 <- t(otu_table(pacu.t1))
OTU.pacu.tf <- t(otu_table(pacu.tf))

###ANCOM pipeline###
##Step 1: set all your parameters
##These should remain the same for all the contrasts we do UNLESS you want to add random effects or covariates. 
##Because we are running contrasts on each species at each time point between treatment, we do no tneed any covariates. 
##Random effects gives us trouble because of singularities (ANCOM uses the nlme package) therefore, all taxa that pull out as
#differentially abundant will also be run through an lmer() from the lme4 package to ensure no difference when Parent.ID is included as a random effect

sample_var <- "Sample.ID"
group_var <- NULL
out_cut <- 0.05  #numerical fraction, below 5% = outlier zeros, above 95% = outlier values
zero_cut <- 0.90 #numerical fraction, taxa with proportion of zeros above 90% are removed
lib_cut <- 1000 #number, removes any samples with less than lib_cut reads
neg_lb <- TRUE

main_var <- "Treatment"
p_adj_method <- "BH" #default: Benjamini-Hochberg procedure
alpha <- 0.05 #level of significance
adj_formula <- NULL
rand_formula <- NULL #add in random effect 

#I just run through the following steps for each of the time points (separate ANCOMs)
#then I combine downstream to make a faceted graph for each species T1 and TF

#M.cap all time points
#feature_table <- OTU.mcap
#meta_data <- sample_data(mcap)

#M cap T1
feature_table <- OTU.mcap.t1
meta_data <- sample_data(mcap.t1)
#M cap TF
feature_table <- OTU.mcap.tf
meta_data <- sample_data(mcap.tf)

levels(meta_data$Treatment) #Ambient first, High Second = High will be positive CLR on graph 
meta_data$Sample.ID <- as.character(meta_data$sample_name.1) #make sure your sample id is chr

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var,
                                     group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM
res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)

#Extract your data frames 

#For all time points together
#resdf.mcap <- as.data.frame(res$out)
#figdf.mcap <- as.data.frame(res$fig$data)

##T1
#Extract your dataframes
resdf.mcap.t1 <- as.data.frame(res$out)
figdf.mcap.t1 <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_t1 = ifelse(is.null(struc_zero), nrow(figdf.mcap.t1), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_t1 = c(0.9 * (n_taxa_t1 -1), 0.8 * (n_taxa_t1 -1), 0.7 * (n_taxa_t1 -1), 0.6 * (n_taxa_t1 -1))
names(cut_off_t1) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.mcap.t1$contrast <- "T1"
figdf.mcap.t1$detected_0.6 <- "240.6"
rownames(figdf.mcap.t1) <- NULL

##Do the same for TF
resdf.mcap.tf <- as.data.frame(res$out)
figdf.mcap.tf <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_tf = ifelse(is.null(struc_zero), nrow(figdf.mcap.tf), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_tf = c(0.9 * (n_taxa_tf -1), 0.8 * (n_taxa_tf -1), 0.7 * (n_taxa_tf -1), 0.6 * (n_taxa_tf -1))
names(cut_off_tf) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.mcap.tf$contrast <- "TF"
figdf.mcap.tf$detected_0.6 <- "89.4"
rownames(figdf.mcap.tf) <- NULL

#Combine your two time points into one dataframe
resdf.mcap <- rbind(resdf.mcap.t1, resdf.mcap.tf)
figdf.mcap <- rbind(figdf.mcap.t1, figdf.mcap.tf)

##Export the data as a data frame
write.csv(resdf.mcap, "ancom/output/ancom_mcap_res.csv")
write.csv(figdf.mcap, "ancom/output/ancom_mcap_figdata.csv")

##Extract the data for use in bubble plots (see figures_ace21.R script)
mcap.true <- filter(resdf.mcap, detected_0.6 == "TRUE")
mcap.otus <- (mcap.true$taxa_id)
mcap.ra <- transform_sample_counts(mcap, function(x) x/ sum(x))
mcap.ancom <- prune_taxa(mcap.otus,mcap.ra)
mcap.ancom.melt <- psmelt(mcap.ancom)
mcap.ancom.melt$top.group <- "differential"
write.csv(mcap.ancom.melt, "ancom/output/mcap_ancom.csv")

##Step 3. Volcano Plot

#Add taxonomic names
tax<-as(tax_table(mcap),"matrix")
#head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
#head(tax)

##Merge taxonomy with your figure data
figdf.mcap <- merge(figdf.mcap,tax, by = "taxa_id")
resdf.mcap <- merge(resdf.mcap, tax, by = "taxa_id")

#Make sure your detected_0.6 level is a numerical vector
figdf.mcap$detected_0.6 <- as.numeric(figdf.mcap$detected_0.6)
#Create annotation data for figure
dat_ann = data.frame(x = min(figdf.mcap$x), y = figdf.mcap$detected_0.6, label = "W[0.6]")

#Order the genera and name only significant genera (all else as OTHER)
x = tapply(figdf.mcap$y, figdf.mcap$Genus, function(x) max(x))
x = sort(x, TRUE)
cut_off #check your cutoff for your intended level - here 0.6
x #look at which taxa are above your cutoff & add to the below code for setting NAs
figdf.mcap$Genus = factor(as.character(figdf.mcap$Genus), levels=names(x))
figdf.mcap$col_genus <- figdf.mcap$Genus

#Set everything that is super low to NA so that we can call them "other"
figdf.mcap$col_genus[figdf.mcap$col_genus != "Alteromonas" & 
                       figdf.mcap$col_genus != "Pseudoalteromonas" &
                       figdf.mcap$col_genus != "Ruegeria" &
                       figdf.mcap$col_genus != "Alteromonadales_unclassified" &
                       figdf.mcap$col_genus != "Reichenbachiella" &
                       figdf.mcap$col_genus != "Rhodobacteraceae_unclassified" &
                       figdf.mcap$col_genus != "uncultured" &
                       figdf.mcap$col_genus != "Prosthecochloris" &
                       figdf.mcap$col_genus != "JGI_0000069-P22_ge"] <- NA


levels(figdf.mcap$col_genus)
# add new factor
figdf.mcap$col_genus <- factor(figdf.mcap$col_genus, levels = c(levels(figdf.mcap$col_genus), "Other"))
# convert NAs to other
figdf.mcap$col_genus[is.na(figdf.mcap$col_genus)] = "Other"

#Make your figure! Realise here that it will plot the two 0.6 level cut-off line
#I remove the second line in illustrator because ggplot is giving me a headache. 
mcap.p <- ggplot(figdf.mcap, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 3) +
  facet_grid(~contrast) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  #geom_hline(yintercept = c(18, linetype = "dashed") + 
  geom_hline(yintercept = figdf.mcap$detected_0.6, linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0.4, color = "orange", parse = TRUE) +
  ggtitle("Montipora capitata")
ggsave("ancom/output/mcap_ancom.pdf")

##All the same steps for Porites compressa! 

##Step 1: set all your parameters
#P comp all
#feature_table <- OTU.pcomp
#meta_data <- sample_data(pcomp)

##If splitting by time point use the following:
##P comp T1
feature_table <- OTU.pcomp.t1
meta_data <- sample_data(pcomp.t1)
##P comp TF
feature_table <- OTU.pcomp.tf
meta_data <- sample_data(pcomp.tf)

levels(meta_data$Treatment) #Ambient first, High Second = High will be positive 
meta_data$Sample.ID <- as.character(meta_data$sample_name.1) #make sure your sample id is chr

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var,
                                     group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM
res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)

##T1
#Extract your dataframes
resdf.pcomp.t1 <- as.data.frame(res$out)
figdf.pcomp.t1 <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_t1 = ifelse(is.null(struc_zero), nrow(figdf.pcomp.t1), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_t1 = c(0.9 * (n_taxa_t1 -1), 0.8 * (n_taxa_t1 -1), 0.7 * (n_taxa_t1 -1), 0.6 * (n_taxa_t1 -1))
names(cut_off_t1) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.pcomp.t1$contrast <- "T1"
figdf.pcomp.t1$detected_0.6 <- "35.4"
rownames(figdf.pcomp.t1) <- NULL

##Do the same for TF
resdf.pcomp.tf <- as.data.frame(res$out)
figdf.pcomp.tf <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_tf = ifelse(is.null(struc_zero), nrow(figdf.pcomp.tf), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_tf = c(0.9 * (n_taxa_tf -1), 0.8 * (n_taxa_tf -1), 0.7 * (n_taxa_tf -1), 0.6 * (n_taxa_tf -1))
names(cut_off_tf) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.pcomp.tf$contrast <- "TF"
figdf.pcomp.tf$detected_0.6 <- "23.4"
rownames(figdf.pcomp.tf) <- NULL

#Combine your two time points into one dataframe
resdf.pcomp <- rbind(resdf.pcomp.t1, resdf.pcomp.tf)
figdf.pcomp <- rbind(figdf.pcomp.t1, figdf.pcomp.tf)

##Export the data as a data frame
write.csv(resdf.pcomp, "ancom/output/ancom_pcomp_res.csv")
write.csv(figdf.pcomp, "ancom/output/ancom_pcomp_figdata.csv")

##Extract the data for use in bubble plots (see figures_ace21.R script)
pcomp.true <- filter(resdf.pcomp, detected_0.6 == "TRUE")
pcomp.otus <- (pcomp.true$taxa_id)
pcomp.ra <- transform_sample_counts(pcomp, function(x) x/ sum(x))
pcomp.ancom <- prune_taxa(pcomp.otus,pcomp.ra)
pcomp.ancom.melt <- psmelt(pcomp.ancom)
pcomp.ancom.melt$top.group <- "differential"
write.csv(pcomp.ancom.melt, "ancom/output/pcomp_ancom.csv")

##Step 3. Volcano Plot

#Add taxonomic names
tax<-as(tax_table(pcomp),"matrix")
#head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
#head(tax)

##Merge taxonomy with your figure data
figdf.pcomp <- merge(figdf.pcomp,tax, by = "taxa_id")
resdf.pcomp <- merge(resdf.pcomp, tax, by = "taxa_id")

#Make sure your detected_0.6 level is a numerical vector
figdf.pcomp$detected_0.6 <- as.numeric(figdf.pcomp$detected_0.6)
#Create annotation data for figure
dat_ann = data.frame(x = min(figdf.pcomp$x), y = figdf.pcomp$detected_0.6, label = "W[0.6]")

#Order the genera and name only significant genera (all else as OTHER)
x = tapply(figdf.pcomp$y, figdf.pcomp$Genus, function(x) max(x))
x = sort(x, TRUE)
cut_off #check your cutoff for your intended level - here 0.6
x #look at which taxa are above your cutoff & add to the below code for setting NAs
figdf.pcomp$Genus = factor(as.character(figdf.pcomp$Genus), levels=names(x))
figdf.pcomp$col_genus <- figdf.pcomp$Genus

#Set everything that is super low to NA so that we can call them "other"
figdf.pcomp$col_genus[figdf.pcomp$col_genus != "Francisellaceae_ge" & 
                        figdf.pcomp$col_genus != "Francisellaceae_unclassified" &
                        figdf.pcomp$col_genus != "Endozoicomonas" &
                        figdf.pcomp$col_genus != "Proteobacteria_unclassified" &
                        figdf.pcomp$col_genus != "Rhodobacteraceae_unclassified" &
                        figdf.pcomp$col_genus != "Bythopirellula" &
                        figdf.pcomp$col_genus != "Sva0996_marine_group"] <- NA


levels(figdf.pcomp$col_genus)
# add new factor
figdf.pcomp$col_genus <- factor(figdf.pcomp$col_genus, levels = c(levels(figdf.pcomp$col_genus), "Other"))
# convert NAs to other
figdf.pcomp$col_genus[is.na(figdf.pcomp$col_genus)] = "Other"

#Make your figure! Realise here that it will plot the two 0.6 level cut-off line
#I remove the second line in illustrator because ggplot is giving me a headache. 
pcomp.p <- ggplot(figdf.pcomp, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 3) +
  facet_grid(~contrast) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  #geom_hline(yintercept = c(18, linetype = "dashed") + 
  geom_hline(yintercept = figdf.pcomp$detected_0.6, linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0.4, color = "orange", parse = TRUE) +
  ggtitle("Porites compressa")
ggsave("ancom/output/pcomp_ancom.pdf")


###Next for Pavona varians
##Step 1: set all your parameters
#P var all
#feature_table <- OTU.pvar
#meta_data <- sample_data(pvar)

#Pvar T1
feature_table <- OTU.pvar.t1
meta_data <- sample_data(pvar.t1)
#Pvar TF
feature_table <- OTU.pvar.tf
meta_data <- sample_data(pvar.tf)

levels(meta_data$Treatment) #Ambient first, High Second = High will be positive 
meta_data$Sample.ID <- as.character(meta_data$sample_name.1) #make sure your sample id is chr

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var,
                                     group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM
res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)

##T1
#Extract your dataframes
resdf.pvar.t1 <- as.data.frame(res$out)
figdf.pvar.t1 <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_t1 = ifelse(is.null(struc_zero), nrow(figdf.pvar.t1), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_t1 = c(0.9 * (n_taxa_t1 -1), 0.8 * (n_taxa_t1 -1), 0.7 * (n_taxa_t1 -1), 0.6 * (n_taxa_t1 -1))
names(cut_off_t1) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.pvar.t1$contrast <- "T1"
figdf.pvar.t1$detected_0.6 <- "128.4"
rownames(figdf.pvar.t1) <- NULL

##Do the same for TF
resdf.pvar.tf <- as.data.frame(res$out)
figdf.pvar.tf <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_tf = ifelse(is.null(struc_zero), nrow(figdf.pvar.tf), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_tf = c(0.9 * (n_taxa_tf -1), 0.8 * (n_taxa_tf -1), 0.7 * (n_taxa_tf -1), 0.6 * (n_taxa_tf -1))
names(cut_off_tf) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.pvar.tf$contrast <- "TF"
figdf.pvar.tf$detected_0.6 <- "93.6"
rownames(figdf.pvar.tf) <- NULL

#Combine your two time points into one dataframe
resdf.pvar <- rbind(resdf.pvar.t1, resdf.pvar.tf)
figdf.pvar <- rbind(figdf.pvar.t1, figdf.pvar.tf)

##Export the data as a data frame
write.csv(resdf.pvar, "ancom/output/ancom_pvar_res.csv")
write.csv(figdf.pvar, "ancom/output/ancom_pvar_figdata.csv")

##Extract the data for use in bubble plots (see figures_ace21.R script)
pvar.true <- filter(resdf.pvar, detected_0.6 == "TRUE")
pvar.otus <- (pvar.true$taxa_id)
pvar.ra <- transform_sample_counts(pvar, function(x) x/ sum(x))
pvar.ancom <- prune_taxa(pvar.otus,pvar.ra)
pvar.ancom.melt <- psmelt(pvar.ancom)
pvar.ancom.melt$top.group <- "differential"
write.csv(pvar.ancom.melt, "ancom/output/pvar_ancom.csv")

##Step 3. Volcano Plot

#Add taxonomic names
tax<-as(tax_table(pvar),"matrix")
#head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
#head(tax)

##Merge taxonomy with your figure data
figdf.pvar <- merge(figdf.pvar,tax, by = "taxa_id")
resdf.pvar <- merge(resdf.pvar, tax, by = "taxa_id")

#Make sure your detected_0.6 level is a numerical vector
figdf.pvar$detected_0.6 <- as.numeric(figdf.pvar$detected_0.6)
#Create annotation data for figure
dat_ann = data.frame(x = min(figdf.pvar$x), y = figdf.pvar$detected_0.6, label = "W[0.6]")

#Order the genera and name only significant genera (all else as OTHER)
x = tapply(figdf.pvar$y, figdf.pvar$Genus, function(x) max(x))
x = sort(x, TRUE)
cut_off #check your cutoff for your intended level - here 0.6
x #look at which taxa are above your cutoff & add to the below code for setting NAs
figdf.pvar$Genus = factor(as.character(figdf.pvar$Genus), levels=names(x))
figdf.pvar$col_genus <- figdf.pvar$Genus

#Set everything that is super low to NA so that we can call them "other"
figdf.pvar$col_genus[figdf.pvar$col_genus != "Bacteria_unclassified"] <- NA


levels(figdf.pvar$col_genus)
# add new factor
figdf.pvar$col_genus <- factor(figdf.pvar$col_genus, levels = c(levels(figdf.pvar$col_genus), "Other"))
# convert NAs to other
figdf.pvar$col_genus[is.na(figdf.pvar$col_genus)] = "Other"

#Make your figure! Realise here that it will plot the two 0.6 level cut-off line
#I remove the second line in illustrator because ggplot is giving me a headache. 
pvar.p <- ggplot(figdf.pvar, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 3) +
  facet_grid(~contrast) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  #geom_hline(yintercept = c(18, linetype = "dashed") + 
  geom_hline(yintercept = figdf.pvar$detected_0.6, linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0.4, color = "orange", parse = TRUE) +
  ggtitle("Pavona varians")
ggsave("ancom/output/pvar_ancom.pdf")

###Finally Pocillopora acuta!
##Step 1: set all your parameters
#P acu all
#feature_table <- OTU.pacu
#meta_data <- sample_data(pacu)

#P acu T1
feature_table <- OTU.pacu.t1
meta_data <- sample_data(pacu.t1)
#P acu TF
feature_table <- OTU.pacu.tf
meta_data <- sample_data(pacu.tf)

levels(meta_data$Treatment) #Ambient first, High Second = High will be positive 
meta_data$Sample.ID <- as.character(meta_data$sample_name.1) #make sure your sample id is chr

pre_pro <- feature_table_pre_process(feature_table, meta_data, sample_var,
                                     group_var,out_cut, zero_cut, lib_cut, neg_lb)
feature_table <- pre_pro$feature_table
meta_data <- pre_pro$meta_data
struc_zero <- pre_pro$structure_zeros

#Step 2 run the ANCOM
res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha,
             adj_formula, rand_formula)

##T1
#Extract your dataframes
resdf.pacu.t1 <- as.data.frame(res$out)
figdf.pacu.t1 <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_t1 = ifelse(is.null(struc_zero), nrow(figdf.pacu.t1), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_t1 = c(0.9 * (n_taxa_t1 -1), 0.8 * (n_taxa_t1 -1), 0.7 * (n_taxa_t1 -1), 0.6 * (n_taxa_t1 -1))
names(cut_off_t1) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.pacu.t1$contrast <- "T1"
figdf.pacu.t1$detected_0.6 <- "202.2"
rownames(figdf.pacu.t1) <- NULL

##Do the same for TF
resdf.pacu.tf <- as.data.frame(res$out)
figdf.pacu.tf <- as.data.frame(res$fig$data)

# Number of taxa except structural zeros
n_taxa_tf = ifelse(is.null(struc_zero), nrow(figdf.pacu.tf), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_tf = c(0.9 * (n_taxa_tf -1), 0.8 * (n_taxa_tf -1), 0.7 * (n_taxa_tf -1), 0.6 * (n_taxa_tf -1))
names(cut_off_tf) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

#Add this info into your figure dataframe for plotting below
figdf.pacu.tf$contrast <- "TF"
figdf.pacu.tf$detected_0.6 <- "225.6"
rownames(figdf.pacu.tf) <- NULL

#Combine your two time points into one dataframe
resdf.pacu <- rbind(resdf.pacu.t1, resdf.pacu.tf)
figdf.pacu <- rbind(figdf.pacu.t1, figdf.pacu.tf)

##Export the data as a data frame
write.csv(resdf.pacu, "ancom/output/ancom_pacu_res.csv")
write.csv(figdf.pacu, "ancom/output/ancom_pacu_figdata.csv")

##Extract the data for use in bubble plots (see figures_ace21.R script)
pacu.true <- filter(resdf.pacu, detected_0.6 == "TRUE")
pacu.otus <- (pacu.true$taxa_id)
pacu.ra <- transform_sample_counts(pacu, function(x) x/ sum(x))
pacu.ancom <- prune_taxa(pacu.otus,pacu.ra)
pacu.ancom.melt <- psmelt(pacu.ancom)
pacu.ancom.melt$top.group <- "differential"
write.csv(pacu.ancom.melt, "ancom/output/pacu_ancom.csv")

##Step 3. Volcano Plot

#Add taxonomic names
tax<-as(tax_table(pacu),"matrix")
#head(tax)
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
#head(tax)

##Merge taxonomy with your figure data
figdf.pacu <- merge(figdf.pacu,tax, by = "taxa_id")
resdf.pacu <- merge(resdf.pacu, tax, by = "taxa_id")

#Make sure your detected_0.6 level is a numerical vector
figdf.pacu$detected_0.6 <- as.numeric(figdf.pacu$detected_0.6)
#Create annotation data for figure
dat_ann = data.frame(x = min(figdf.pacu$x), y = figdf.pacu$detected_0.6, label = "W[0.6]")

#Order the genera and name only significant genera (all else as OTHER)
x = tapply(figdf.pacu$y, figdf.pacu$Genus, function(x) max(x))
x = sort(x, TRUE)
cut_off #check your cutoff for your intended level - here 0.6
x #look at which taxa are above your cutoff & add to the below code for setting NAs
figdf.pacu$Genus = factor(as.character(figdf.pacu$Genus), levels=names(x))
figdf.pacu$col_genus <- figdf.pacu$Genus

#Set everything that is super low to NA so that we can call them "other"
figdf.pacu$col_genus[figdf.pacu$col_genus != "Pseudoalteromonas" &
                       figdf.pacu$col_genus != "Synechococcus_CC9902" & 
                       figdf.pacu$col_genus != " Filomicrobium"] <- NA


levels(figdf.pacu$col_genus)
# add new factor
figdf.pacu$col_genus <- factor(figdf.pacu$col_genus, levels = c(levels(figdf.pacu$col_genus), "Other"))
# convert NAs to other
figdf.pacu$col_genus[is.na(figdf.pacu$col_genus)] = "Other"

#Make your figure! Realise here that it will plot the two 0.6 level cut-off line
#I remove the second line in illustrator because ggplot is giving me a headache. 
pacu.p <- ggplot(figdf.pacu, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 3) +
  facet_grid(~contrast) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  #geom_hline(yintercept = c(18, linetype = "dashed") + 
  geom_hline(yintercept = figdf.pacu$detected_0.6, linetype = "dashed") +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0.4, color = "orange", parse = TRUE) +
  ggtitle("Pocillopora acuta")
ggsave("ancom/output/pacu_ancom.pdf")

############# Algal phylotype??
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
#subset to only the experimental corals (no field)
bac.exp <- subset_samples(Bac.seq, Type == "sample")
##Fix sample data so that all T0 are labelled as "Ambient"
#Concatenate & make a new data column to combine all data at T0 as "ambient"
sample_data(bac.exp)['Treat.Time'] <- paste(sample_data(bac.exp)$Treatment, sample_data(bac.exp)$Time.Point)
View(as(sample_data(bac.exp), "data.frame"))
#Revalue
sample_data(bac.exp)$Treat.Time <- as.factor(sample_data(bac.exp)$Treat.Time)
levels(sample_data(bac.exp)$Treat.Time)
sample_data(bac.exp)$Treat.Time <- revalue(sample_data(bac.exp)$Treat.Time, c("Ambient T0" = "Ambient", "Ambient T1" = "Ambient", "Ambient TF" = "Ambient", "High T0" = "Ambient", 
                                                                              "High T1" = "High", "High TF" = "High"))
levels(sample_data(bac.exp)$Treat.Time)

mcap.exp <- subset_samples(bac.exp, Species == "Montipora_capitata")
mcap.exp.ra <- transform_sample_counts(mcap.exp, function(x) x/ sum(x))
bac.exp.ra <- transform_sample_counts(bac.exp, function(x) x/ sum(x))
endo <- subset_taxa(bac.exp.ra, Genus == "Endozoicomonas")
merged.endo <- tax_glom(endo, "Genus", NArm = FALSE)
endo.melt <- psmelt(merged.endo)
sum.endo <- ddply(endo.melt, c("Species", "Time.Point", "Treat.Time", "Clade"), summarise,
                     N = length(Abundance), 
                     mean = mean(Abundance),
                     sd = sd(Abundance), 
                     se = sd/sqrt(N)
)
sum.endo <- ddply(endo.melt, c("Species", "Time.Point", "Treat.Time"), summarise,
                 N = length(Abundance), 
                 mean = mean(Abundance),
                 sd = sd(Abundance), 
                 se = sd/sqrt(N)
)
ggplot(sum.endo, aes(x = Time.Point, y = mean, fill = Treat.Time)) +
  geom_point(aes(color = Treat.Time), size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = Treat.Time), width = 0.1) +
  scale_color_manual(name = "Treat.Time", values = c("#0842E2", "#DE0505")) +
  facet_grid(cols = vars(Species)) +
  theme_bw()
ggsave("ancom/output/allspecies_treat_endo.pdf")

bac.exp.ra <- transform_sample_counts(bac.exp, function(x) x/ sum(x))
pseudo <- subset_taxa(bac.exp.ra, Genus == "Pseudoalteromonas")
merged.pseudo <- tax_glom(pseudo, "Genus", NArm = FALSE)
pseudo.melt <- psmelt(merged.pseudo)
sum.pseudo <- ddply(pseudo.melt, c("Species", "Time.Point", "Treat.Time", "Clade"), summarise,
                  N = length(Abundance), 
                  mean = mean(Abundance),
                  sd = sd(Abundance), 
                  se = sd/sqrt(N)
)
sum.pseudo <- ddply(pseudo.melt, c("Species", "Time.Point", "Treat.Time"), summarise,
                  N = length(Abundance), 
                  mean = mean(Abundance),
                  sd = sd(Abundance), 
                  se = sd/sqrt(N)
)
ggplot(sum.pseudo, aes(x = Time.Point, y = mean, fill = Treat.Time)) +
  geom_point(aes(color = Treat.Time), size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = Treat.Time), width = 0.1) +
  scale_color_manual(name = "Treat.Time", values = c("#0842E2", "#DE0505")) +
  facet_grid(cols = vars(Species)) +
  ggtitle("Pseudoalteromonas") +
  theme_bw()
ggsave("ancom/output/allspecies_treat_endo.pdf")


##Pavona -  C1 vs. C27 
##What about pavona vs acuta 


