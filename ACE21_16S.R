# load in 16S to phyloseq

library(phyloseq)
library(readr)



# this sorft of worked https://www.cmich.edu/colleges/se/biology/fac_staff/Learman%20Research%20Lab/ANT_16S_R.txt
# #load in tax and otu
# biom1<-read.csv("Data/TT17_3000_summer_16S-pipeline_outputs/Results/main/raw_abundanceTable_100.shared")
# mp0<-import_biom(biom1)
# tax_table(mp0) <- tax_table(mp0)[, 1:6] #keeping only 6 tax levels apparently this fixes a bug


#read in sample data
MetaData16<-read.csv("Data/Meta_16s_new.csv") #run in physio script first
MetaData16$Parent.ID<-as.factor(as.character(MetaData16$Parent.ID))
MetaData16$Tank.num<-as.factor(as.character(MetaData16$Tank.num))

sam0 <- MetaData16
sam1 <- as.matrix(sam0[, -1])
rownames(sam1) <- sam0$sample_name
sam <- sample_data(data.frame(sam1))


#load in OUT:raw_abundanceTable_100.shared

OTU2<-read.table("Data/TT17_3000_summer_16S-pipeline_outputs/raw_abundanceTable_100_shared.csv", sep=',', header=T)

#from its
#otu
otu5 <- as.matrix(OTU2[, -1])
rownames(otu5) <- OTU2$sample_name
otu <- otu_table(otu5, taxa_are_rows = FALSE)


#tax
TAX<- read.csv("Data/TT17_3000_summer_16S-pipeline_outputs/raw_consensusClassification_100_taxonomy.csv", colClasses = "character")
tax1 <- as.matrix(TAX[, -1], dimnames = list(TAX$OTU, colnames(TAX[-1])))
rownames(tax1) <- tax0$OTU
tax <- tax_table(tax1)




# Read the data into phyloseq
Bac.seqAll = phyloseq(otu, tax,sam) #THIS WORKS
Bac.seqAll

save(Bac.seqAll, file = "Data/Bac.seqAll_phyloseq.RData")

## check out data
ntaxa(Bac.seqAll)  #num taxa
nsamples(Bac.seqAll)   #num samples
sample_names(Bac.seqAll)[1:300] #samp names
 #sampsy<-sample_names(physeq)[1:300]
#write.csv(sampsy,"sampsy.csv")

# OK Now you need to do the steps in the pipeline.
#1 filter out taxa: get rid of unknown, mitochondria and chloroplast
Bac.seqAll<-subset_taxa(Bac.seqAll, Family !="Mitochondria")
Bac.seqAll<-subset_taxa(Bac.seqAll, Order !="Chloroplast")
Bac.seqAll<-subset_taxa(Bac.seqAll, Taxonomy !="unknown")
Bac.seqAll<-subset_taxa(Bac.seqAll, Taxonomy !="Archaea")

# now get rid of things not represnted 2x or more
#Prune singletons (these are reads that are only found once)
Bac.seqAll.prune <- prune_taxa(taxa_sums(Bac.seqAll) > 2, Bac.seqAll)

  save(Bac.seqAll, file = "Data/Bac.seqAll.prune.RData")
  












####################
#not used here but great tutorial with images https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

#load in 3000 rare instead
#read in sample data
MetaData16<-read.csv("Data/Meta16s_3k.csv") #run in physio script first
MetaData16$Parent.ID<-as.factor(as.character(MetaData16$Parent.ID))
MetaData16$Tank.num<-as.factor(as.character(MetaData16$Tank.num))

sam0 <- MetaData16
sam1 <- as.matrix(sam0[, -1])
rownames(sam1) <- sam0$sample_name
sam <- sample_data(data.frame(sam1))

#load in OUT:raw_abundanceTable_100.shared

OTU3k<-read.table("Data/TT17_3000_summer_16S-pipeline_outputs/abundance_table_100_shared.csv", sep=',', header=T)

#from its  
#otu
otu3k1 <- as.matrix(OTU3k[, -1])
rownames(otu3k1) <- OTU3k$sample_name
otu <- otu_table(otu3k1, taxa_are_rows = FALSE)


#tax table annotations_100.taxonomy.csv
TAX<- read.csv("Data/TT17_3000_summer_16S-pipeline_outputs/annotations_100.taxonomy.csv", colClasses = "character")
tax1 <- as.matrix(TAX[, -1], dimnames = list(TAX$OTU, colnames(TAX[-1])))
rownames(tax1) <- TAX$OTU
tax <- tax_table(tax1)


# Read the data into phyloseq
Bac.seq2 = phyloseq(otu, tax,sam) #THIS WORKS
Bac.seq
Bac.seq.df <- sample_data(Bac.seq)

save(Bac.seq, file = "Data/Bac.seq_phyloseq.RData")


## check out data
ntaxa(Bac.seq)  #num taxa
nsamples(Bac.seq)   #num samples
sample_names(Bac.seq)[1:300] #samp names
 #sampsy<-sample_names(physeq)[1:350]
#write.csv(sampsy,"sampsy.csv")




