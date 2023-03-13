## Post dada2 Data Analysis
## Initial exploration of data (PCoA)
## 16s

## Set directory & check files
setwd("/Users/mjp23/Desktop/RootAgePostDada2")

list.files("/Users/mjp23/Desktop/RootAgePostDada2")

### Load required packages ###
library(ade4)
library(RColorBrewer)
library(phyloseq)
library(plyr)
library(gdata)
library(ggplot2)
library(vegan)
library(agricolae)
library(tidyverse)
library(gtools)

##install.packages("janitor")
library(janitor)

## A) Import and organize the data

## Import data
tax1 <- read.table(file="Ph_RootAge16s_ESV_Taxonomy.txt", header=T, row.names=1)
tax1 ## looks good
tax <- tax1[-c(8)] ## remove sequence column

asv_tab <- read.table(file="Ph_RootAge16s_ESV_Abund_Table.txt", header=T, row.names=1) 
asv_tab ## asv don't need to transpose, but need to remove blanks (95,96,112) and 2nd order roots
data.frame(rownames(asv_tab)) ## want to remove rows that are blanks in metadata file and 2nd order roots
asv <- asv_tab[-c(12,40,48,53, 59,77,95,96,100,112),] 
asv ## looks good; no blank samples and no 2nd order roots

meta <- read.csv("Ph_RootAge16s_metadata.csv", header=T, row.names=1)
meta ## need to remove last row of blank
view(meta)
metadat <- meta[-c(12,40,48,53,59,77,95,96,100,112),]
metadat ## looks good; blanks are removed
metadat$Root.Age <- as.numeric(as.character(metadat$Root.Age))
is.numeric(metadat$Root.Age) 

## Create new columns in metadata
metadat$Vine <- paste(metadat$Row, metadat$Block) ## created a new column, "Vine"
metadat$YoungvsOld <- ifelse(metadat$Root.Age<11,"Young", "Old") ## Creates a new, old vs. young based on Volder et al.
metadat$YoungvsOld

metadat ## looks good

## Remove problematic taxa
taxon<- tax
taxon[is.na(taxon)]<-"Unclassified"
taxon ## Now has "Unclassified" instead of NA

euk <- subset(taxon, Kingdom == "Eukaryota")
unclass <- subset(taxon, Phylum == "Unclassified")
mito <- subset(taxon, Family == "Mitochondria")
arch <- subset(taxon, Kingdom == "Archaea")
chloro <- subset(taxon, Class == "Chloroplast")

rownames(euk)
rownames(unclass)
rownames(mito)
rownames(arch)
rownames(chloro)

un.1<-union(rownames(euk),rownames(unclass))
un.2<-union(rownames(mito),un.1)
un.3<-union(rownames(arch),un.2)
un.4 <-union(rownames(chloro), un.3)
un.4

red.taxon <- taxon[-which(rownames(taxon) %in% un.4),]
red.taxon


## Match asv table with new taxonomy file
asv.filt<-asv[,which(colnames(asv) %in% rownames(red.taxon))]
asv.filt ## just looking at what comes up

## Calculate reads per sample
rowSums(asv.filt) ## 2, 23, 31, 32, 89, 105, 106, 107 have no samples

## Determine minimum number of reads in your data set
min(rowSums(asv.filt))

## Cut samples with too few reads
n2 <- names(which(rowSums(asv.filt) > 1000)) 
asv.re<-asv.filt[which(rownames(asv.filt) %in% n2),]
min(rowSums(asv.re)) ## min now 1377 reads

metadat.red <-metadat[which(rownames(metadat) %in% rownames(asv.re)),] ## match metadata with removed samples

## Rarefy to obtain even numbers of reads by sample 
set.seed(336)
asv.r<-rrarefy(asv.re, 1355)
asv.r

rowSums(asv.r)

## Convert asvs to percentages
asv.perc<-asv.r/rowSums(asv.r)*100
asv.perc

## Set factors as factors before analysis
# For whole dataset
metadat.red$Block <-as.factor(metadat.red$Block)
metadat.red$Row <- as.factor(metadat.red$Row)
metadat.red$Root.type <- as.factor(metadat.red$Root.type)
metadat.red$Age.testing.or.variation.testing. <- as.factor(metadat.red$Age.testing.or.variation.testing.)
metadat.red$Vine <- as.factor(metadat.red$Vine)
metadat.red$Mass..mg.<- as.factor(metadat.red$Mass..mg.)
metadat.red$Color <-as.factor(metadat.red$Color)
metadat.red$AgeCat <- as.factor(metadat.red$AgeCat)
metadat.red$YoungvsOld<-as.factor(metadat.red$YoungvsOld)


## Counting root and vine distribution of full data set
metadat.red %>% 
  group_by(metadat.red$Vine) %>% 
  summarize(number_rows=n())
metadat.red %>%
  group_by(metadat.red$Root.Age) %>%
  summarize(number_rows=n())

tabyl(metadat.red, Root.Age, Vine) ## table shows root age, vine, and number of samples distribution

## CREATE PHYLOSEQ OBJECT
## Group the asv, metadata, and taxonomy files
asv.perc<-as.data.frame(asv.perc)
asv.perc

As.3.ASV <- otu_table(asv.perc, taxa_are_rows=F)
As.3.Meta <- sample_data(metadat.red,errorIfNULL=TRUE)
As.3.Taxo <- tax_table(as.matrix(red.taxon), errorIfNULL=TRUE)
As.3.Phylo <- phyloseq(As.3.ASV, As.3.Meta,As.3.Taxo)

nsamples(As.3.Phylo)

## PCOA Plots 
Rootage16s.ord <- ordinate(As.3.Phylo, "PCoA", "bray")

p1 = plot_ordination(As.3.Phylo, Rootage16s.ord, type="taxa", color="Phylum", title="taxa")
print(p1)

p2= plot_ordination(As.3.Phylo, Rootage16s.ord, type="samples", color="Vine", shape="YoungvsOld") + geom_point(size=4) +
labs(title = " PCoA 16s: All sampled roots (n = 91)") +
  theme_bw ()
print(p2) ## PCOA including information about vine

## Data shows separation of roots by color (vine); therefore, subsequent analysis comparing young vs. old roots will 
## factor in the variation caused by vine; furthermore, vine 21A only has one root, and it will be removed.