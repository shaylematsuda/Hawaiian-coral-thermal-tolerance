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

#Tell phyloseq what is what otu and tax tables need to be matrix
OTU = otu_table(OTU3, taxa_are_rows = T)
TAX = tax_table(TAX2)

# Read the data into phyloseq
Bac.seq = phyloseq(otu, tax,sam) #THIS WORKS
Bac.seq

save(Bac.seq, file = "Data/Raw_bacteria_phyloseq.RData")


## check out data
ntaxa(physeq)  #num taxa
nsamples(physeq)   #num samples
sample_names(physeq)[1:300] #samp names
# sampsy<-sample_names(physeq)[1:300]
#write.csv(sampsy,"sampsy.csv")
