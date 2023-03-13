### Meredith Persico
### Root Age 16S Pre-Processing
## Started Analysis 10/13/21 in ROAR; R graphical interface
## [1] "/storage/home/m/mxp1191/R/x86_64-pc-linux-gnu-library/4.1"   


## Set the working directory 

list.files() #This command shows files in your working directory
getwd()
setwd("/storage/work/mxp1191/Rootage_16S")
list.files()

###install packages ###
install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("tibble")
install.packages("ggfortify")
install.packages("dplyr")
install.packages("Biostrings")
install.packages("ape")
install.packages("tidyr")

# Load packages
library("dada2")
library("ShortRead")
library("tibble")
library("phyloseq")
library("ggplot2")
library("vegan")
library("dplyr")
library("Biostrings")
library("ShortRead")
library("ape")
library("tidyr")


#### ASSIGN PATH ####

# CHANGE PATH to the directory containing the fastq files after unzipping.
path <- ("/storage/work/mxp1191/Rootage_16S")
path

list.files(path)

#### READ IN FILES ####

# First we read in the names of the fastq files, and perform some string manipulation to get lists of the 
# forward and reverse fastq files in matched order

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
#last number is the set location between _
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 5)

sample.names

#### EXAMINE QUALITY PROFILES OF FORWARD AND REVERSE READS #### 

## forward reads ##
# Poor but not terrible quality- trim at 220

plotQualityProfile(fnFs[3:6])
plotQualityProfile(fnFs[8:11])
plotQualityProfile(fnFs[13:16])

## reverse reads ##
# Poor but not terrible quality- trim at 180

plotQualityProfile(fnRs[3:6])
plotQualityProfile(fnRs[8:11])
plotQualityProfile(fnRs[13:16])


#### FILTERING AND TRIMMING #### 

# Assign the filenames for the filtered fastq.gz files

# Place filtered files in filtered/ subdirectory

filt_path <- file.path(path, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

filtFs
filtRs

# Filter forward and reverse reads
# On Windows set multithread=FALSE
# These are standard filtering parameters
# maxEE sets the max expected error - better than avg of quality scores
# For common ITS amplicon strategies, it is undesirable to truncate reads to a fixed length due to the large amount of length variation at that locus. 
#That is OK, just leave out truncLen. 
#Make sure you removed the forward and reverse primers from both the forward and reverse reads though!

#NOTE unlike ITS can trim

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
out

#### LEARN THE ERROR RATES #### 

## The program has to estimate the sequencing error rate, to try and 
## distinguish errors from true Exact Sequence Variants
# Iterative process of comparing possible error structures
# STEP TAKES LONG TIME

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

## Visualize errors
# plots each possible error transition between bases
# points are observed error rates, black line is estimated error, red line is expected error

plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)


#### DEREPLICATION (equivalent of unique.seqs) #### 

# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#### SAMPLE INFERENCE #### 

# Infer the sequence variants in each sample #

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspecting the dada-class object returned by dada:#
# In other words, how many real sequence variants do you have? #

dadaFs[[1]]
dadaRs[[1]]

#### MERGE PAIRED READS #### 

# Merge the denoised forward and reverse reads: #

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#### CONSTRUCT SEQUENCE TABLE #### 

#We can now construct a sequence table, a higher-resolution version of the OTU table produced by traditional methods.#

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

#### REMOVE CHIMERAS #### 

# Remove chimeric sequences:

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Percent of sequences that are non-chimeric
## This differs from number of sequence variants removed, which may be high
sum(seqtab.nochim)/sum(seqtab)

## 0.608228

#### TRACK READS THROUGH PIPELINE #### 

getN <- function(x) sum(getUniques(x))

# Combine sequence reductions by sample from each step #
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

# Show sequence reductions throughout pipeline
track

# Save tracked files for record
write.csv(track, file="tracked_seq_Rootage_16S")

#### ASSIGN TAXONOMY #### 

# Make sure to include PATH info for Silva database
taxa <- assignTaxonomy(seqtab.nochim, "/storage/work/mxp1191/Rootage_16S/silva_nr_v128_train_set.fa", multithread=TRUE)
head(taxa)

## Optionally we can add species assignments due to higher resolution ##

taxa_sp <- addSpecies(taxa, "/storage/work/mxp1191/Rootage_16S/silva_species_assignment_v128.fa")
taxa <- taxa_sp ## Originally created "taxa_sp" but don't want a separate one with just species
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save the important end results for easy upload later
save(taxa, seqtab.nochim, file = "Rootage_16s_taxa_seqtab.Rdata")

#################################################################

####Make PS Object####
library(phyloseq)
library(dplyr)
library(tidyverse)

# three files for phyloseq, esv table (otu table), taxonomy table, metadata
#import the metadata file
Map<-as.data.frame(read.csv(file="/storage/work/mxp1191/Rootage_16S/RootAge_Metadata.csv", header= TRUE, sep=",",row.names=1))
Map

# making phyloseq requires common linker to merge the files into one object
# take a look at each file individual
seqtab.nochim # you'll notice that the column names are sequences, not ESVs
taxa # you'll notice that the row names are sequences
Map # row names are sample names

##rownames(seqtab.nochim)==c(270:372) #double check all rownames in order

##rownames(seqtab.nochim)=c(1:103)#replace with sampleID

# we need to change the sequences to ESVs  for both the seqtab.nochim (esv table) and taxa table
# A. Duplicate your final ESV table #
#seqtab.nochim.mod<-seqtab.nochim.filt
seqtab.nochim.mod<-as.data.frame(seqtab.nochim)

## B. Replace column names for ESV, table with ESV numbers 
ncol(seqtab.nochim.mod)
# 19610 ESVs
fin.ESV.names<-paste("B_asv",sep = "",rep(1:ncol(seqtab.nochim.mod)))

colnames(seqtab.nochim.mod)<-fin.ESV.names

# open up seqtab.nochim.mod, sequences are now replaced by ESVs (easier to manipulate later)
seqtab.nochim.mod ## looks good now

esv.table<-otu_table(seqtab.nochim.mod, taxa_are_rows=FALSE)
# this formats the esv table into an "otu_table", a requirement for phyloseq objects
esv.table

## C. Add ESV names to Taxonomy table and then write to file
# We will add the ESVs from the ESV table to the taxonomy table 
# then move the column of ESVs into the rows

fin.ESV.names<-colnames(esv.table) # puts the ESVs into a character string

taxa.mod<-cbind(taxa, fin.ESV.names) # adds the character string of ESVs to the taxa table
# open up taxa.mod

taxa.df<-rownames_to_column(as.data.frame(taxa.mod), var="sequence") 
# puts sequences into a column, because we want to keep them. Column name is sequence

taxa.evs.final<-column_to_rownames(taxa.df, var="fin.ESV.names")
# this moves the column "fin.ESV.names" to the rows
taxa.evs.final

# now we have to format the taxa table into a tax_table
tax.evs.final.table<-tax_table(taxa.evs.final)
# we get an error saying ... coerce to matrix manually

tax.evs.final.mat<-as.matrix(taxa.evs.final)
tax.evs.final.table<-tax_table(tax.evs.final.mat)
head(taxa_names(tax.evs.final.table))

##D. Now we will build the phyloseq object. To do so, mapping file must be changed
length(sample_names(esv.table))

Map2=sample_data(Map)

sample_names(Map2)
sample_names(esv.table)
#NO SAMPLES MISSING 

sample_names(Map2)<- sample_names(esv.table)

ps.16s.ALL <- phyloseq(esv.table, tax.evs.final.table, Map2)

saveRDS(ps.16s.ALL, file="Raw.16s.all.rds")

#Now, let's make sure that the correct information is included in our phyloseq object.
#If didn't load in, would just give you "NULL"


ps.16s.ALL@otu_table #Should include ESV info
ps.16s.ALL@tax_table #Should include taxonomic info
ps.16s.ALL@sam_data #Should include sample(MAP) info

#Export ps object for downstream upload and processing

save(ps.16s.ALL, file = "PhyloSeqRootAge16s.Rdata")

#Now export the taxonomy table for records
write.table(tax.evs.final.table, file = "Ph_RootAge16s_ESV_Taxonomy.txt")
## Downloaded this file and moved the sequences to last column

# Write ESV table to file
write.table(seqtab.nochim.mod, file="Ph_RootAge16s_ESV_Abund_Table.txt")
taxtab <- read.table(file="Ph_RootAge16s_ESV_Abund_Table.txt")

#Write Sample data to file (transformed)
write.table(Map2,file="Ph_RootAge16s_metadata.txt")

######################################################################################
#### ALT WRITE OUTPUTS TO STOREABLE FILES ####

### Write to CSVs -- This is what I tried first, but then realized that the format was not correct 10/27/21
##write.csv(taxa, file="TAXA_RootAge16s.csv")
##write.csv(seqtab.nochim, file="seqtab.nochim.RootAge16s.csv")

## FROM DADA2 MICROBIAL ANALYSIS LAB ##

# Duplicate your final ESV table #
seqtab.nochim.mod<- seqtab.nochim

## Replace column names for ESV table with ESV numbers 
ncol(seqtab.nochim.mod)
fin.ESV.names<- paste("ESV",sep = "",rep(1:ncol(seqtab.nochim.mod)))

colnames(seqtab.nochim.mod)<-fin.ESV.names

# Write ESV table to file
write.table(seqtab.nochim.mod, file="RootAge16s_ESV_Abund_Table.txt")

## Add ESV names to Taxonomy table and then write to file ##

taxa.mod<-cbind(fin.ESV.names,taxa)

write.table(taxa.mod, file = "RootAge16s_ESV_Taxonomy.txt")


### More with PhyloSeq

####Plot of %bacteria####
ps.perc<-transform_sample_counts(ps.16s.ALL, function(x) x / sum(x)*100) #This calculates the percentage of each read in each sample

head(as.data.frame(ps.perc@tax_table))

Euk <- subset_taxa(ps.perc, Kingdom == "Eukaryota")
Euk.sum<-rowSums(otu_table(Euk))
Arch <- subset_taxa(ps.perc, Kingdom == "Archaea")
Arch.sum<-rowSums(otu_table(Arch))
Chloro <- subset_taxa(ps.perc, Class == "Chloroplast")
Chloro.sum<-rowSums(otu_table(Chloro))
Mito <- subset_taxa(ps.perc, Family == "Mitochondria")
Mito.sum<-rowSums(otu_table(Mito))


###All left are classified as bacteria###
bacter<-100-(Euk.sum+Arch.sum+Chloro.sum+Mito.sum)

met<-as.data.frame(c(1:103))

E<-cbind(met,Euk.sum) 
E$type<-"1Eukaryota"
colnames(E)<-c("ID","Value","Type")

A<-cbind(met,Arch.sum) 
A$type<-"3Archaea"
colnames(A)<-c("ID","Value","Type")

C<-cbind(met,Chloro.sum) 
C$type<-"4Chloroplasts"
colnames(C)<-c("ID","Value","Type")

M<-cbind(met,Mito.sum) 
M$type<-"2Mitochondria"
colnames(M)<-c("ID","Value","Type")

B<-cbind(met,bacter) 
B$type<-"5Bacteria"
colnames(B)<-c("ID","Value","Type")

tab<-rbind(B,E,A,C,M)

min(B)

quartz()
ggplot(tab, aes(fill=Type, y=Value, x=ID)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=c("grey10","grey20","grey40","grey60","black"))

write.csv(tab,"percentbacteria.csv")


#### create bacteria only PS object ####

taxtab<-as.data.frame(ps.16s.ALL@tax_table)

taxtab[is.na(taxtab)]<-"unclassified"

tax.tab.fin<-as.matrix(taxtab)
ps.16s.ALL@tax_table<-tax_table(tax.tab.fin)

filt1 <- prune_taxa(taxtab$Kingdom != c("Eukaryota","Archaea"), ps.16s.ALL)
taxtab2<-as.data.frame(filt1@tax_table)

filt2 <-prune_taxa(taxtab2$Class != "Chloroplast", filt1)
taxtab3<-as.data.frame(filt2@tax_table)


ps.16s.bac <-prune_taxa(taxtab3$Family != "Mitochondria", filt2)
#from 17600 to 17025 taxa, or 3.3% of taxa

save(ps.16s.bac,file = "PSobj16sBacOnly.Rdata")

