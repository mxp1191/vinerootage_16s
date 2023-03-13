### Meredith Persico: mxp1191@psu.edu
### Code for: "The age of absorptive roots impacts root-adjacent microbial composition in grapevines"
### This code covers differential abundance analysis using ancombc for 16s data


## Set directory & check files
setwd("/Users/mjp23/Desktop/RootAgePostDada2")

list.files("/Users/mjp23/Desktop/RootAgePostDada2") ## insert your directory

### Load required packages ###
### Install packages as needed ###
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
library(ANCOMBC)
library(microbiome)
library(pheatmap)

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
asv <- asv_tab[-c(12,33,40,48,53, 59,77,95,96,100,112),] ### removed sample 16s33 on 3/16/22 (vine with one root)
asv ## looks good; no blank samples, other removed samples above (48,77,12,100,59,40,53 are second order roots)

meta <- read.csv("Ph_RootAge16s_metadata.csv", header=T, row.names=1)
meta ## need to remove last row of blank
metadat <- meta[-c(12,33,40,48,53,59,77,95,96,100,112),] ## removed sample 16s33 on 3/16/22 (vine with one root)
metadat ## looks good; blanks are removed
metadat$Root.Age <- as.numeric(as.character(metadat$Root.Age)) ## changes root age to numeric value
is.numeric(metadat$Root.Age) ## confirms yes root age is numeric

## Create new columns in metadata
metadat$Vine <- paste(metadat$Row, metadat$Block) ## created a new column, "Vine"
metadat$AgeCat <- ifelse(metadat$Root.Age<13,"Young", "MiddleAge") ## Create new column for age category ("young" and "middle age")
metadat$AgeCat <- ifelse(metadat$Root.Age<27, metadat$AgeCat, "Old") ## Create new column for age category ("old")
metadat$AgeCat ## check categories are correct
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

test <- subset(red.taxon, Kingdom == "Eukaryota") ## test to make sure they were removed
test2 <- subset(red.taxon, Phylum == "Unclassified") ## test to make sure they were removed

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
min(rowSums(asv.re)) ## min now 1355 reads; 11 sample difference from asv.filt

metadat.red <-metadat[which(rownames(metadat) %in% rownames(asv.re)),] ## match metadata with removed samples


## Split up age and variation roots from metadata and asv.perc files
ancom.asv.fac <- cbind(metadat.red, asv.re)
ancom.asv.age <- subset(ancom.asv.fac, Age.testing.or.variation.testing.%in% c("Age")) ## split into just age-related roots
ancom.asv.var <- subset(ancom.asv.fac, Age.testing.or.variation.testing.%in% c("Variation")) ## split to just variation-related roots

## "age" data set; 44 observations
ancom.age.met <- ancom.asv.age[,1:11] ## new metadata
ancom.asv.age.noperc <- ancom.asv.age[,12:18039] ## new ASV

## removing the residual zero asvs from when age and var datasets were combined
test <- t(ancom.asv.age.noperc)
rowSums(test)
no_zeroes <- test[rowSums(test)>0,] 
no.zero.taxa <- red.taxon[which(rownames(red.taxon) %in% rownames(no_zeroes)),]
asv.nz <- t(no_zeroes) ## asv
meta.nzes <- metadat.red[which(rownames(metadat.red) %in% rownames(asv.nz)),]


## variation data set
ancom.asv.var.met <- ancom.asv.var[,1:12] ## new metadata
ancom.asv.var.noperc <- ancom.asv.var[,13:18040] ## new ASV

testvar <- t(ancom.asv.var.noperc)
rowSums(testvar)


ancom.age.met$Vine <- as.factor(ancom.age.met$Vine)
ancom.age.met$YoungvsOld <-as.factor(ancom.age.met$YoungvsOld)

#### REMOVED SAMPLES WITH LOW READS <1000 BUT DID NOT RAREFY ####

## create phyloseq object of non-normalized data; includes taxa from variation and "just age" root datasets
ancomb.16s.ASV.age <- otu_table(ancom.asv.age.noperc, taxa_are_rows=F)
ancomb.16s.Meta.age <- sample_data(ancom.age.met,errorIfNULL=TRUE)
ancomb.16s.Taxo.age <- tax_table(as.matrix(red.taxon), errorIfNULL=TRUE)
ancomb.16s.Phylo.age <- phyloseq(ancomb.16s.ASV.age, ancomb.16s.Meta.age,ancomb.16s.Taxo.age) 
ancomb.16s.Phylo.age

save(ancomb.16s.Phylo.age,
     file = "ancom.16s.Phylo.age.Rdata")

#### phyloseq object of non-normalized data and no residual taxa from the "variation" roots
### use this phyloseq object for ancombc, because does not contain extraneous taxa 
nz.16s.ASV.age <- otu_table(asv.nz, taxa_are_rows=F)
nz.16s.Meta.age <- sample_data(meta.nzes,errorIfNULL=TRUE)
nz.16s.Taxo.age <- tax_table(as.matrix(no.zero.taxa, errorIfNULL=TRUE))
nz.16s.Phylo.age <- phyloseq(nz.16s.ASV.age, nz.16s.Meta.age, nz.16s.Taxo.age) 
nz.16s.Phylo.age


########### ANCOMBC: ASV LEVEL

### Step 1: zero_cut threshold is ".99" which means it includes all samples
### and even taxa that are rare.

out = ancombc(phyloseq = nz.16s.Phylo.age, formula = "YoungvsOld", 
              p_adj_method = "holm", zero_cut = 0.99, lib_cut = 1000, 
              group = "YoungvsOld", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

res = out$res
res_global = out$res_global
tab_p = res$p_val ## p value
tab_q = res$q ## adjusted p values


tab_diff = res$diff_abn
as.data.frame(tab_diff)
count(tab_diff, tab_diff$YoungvsOld == "TRUE")
true <- subset(tab_diff, tab_diff$YoungvsOldYoung == "TRUE")
as.data.frame(true)
rownames(tab_diff) 

acomb.true.tax<-red.taxon[which(rownames(red.taxon) %in% rownames(true)), 2:7]
as.data.frame(acomb.true.tax)
tab_w = res$W

out_feat <- out$feature_table
out$samp_frac
out$zero_ind
out$samp_frac
out$resid
out$delta_wls

### Step 2: zero_cut threshold is now "0.75" which means it includes only
### taxa appearing in > 25% of samples

out.40 = ancombc(phyloseq = nz.16s.Phylo.age, formula = "YoungvsOld", 
                 p_adj_method = "holm", zero_cut = 0.75, lib_cut = 1000, 
                 group = "YoungvsOld", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                 max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

res.40 = out.40$res
tab_p40 = res40$p_val ## p value
tab_q40 = res40$q ## adjusted p values


tab_diff40 = res.40$diff_abn
as.data.frame(tab_diff40)
count(tab_diff40, tab_diff40$YoungvsOld == "TRUE")
true40 <- subset(tab_diff40, tab_diff40$YoungvsOldYoung == "TRUE")
as.data.frame(true40)
rownames(tab_diff40) 

acomb.true.tax40<-red.taxon[which(rownames(red.taxon) %in% rownames(true40)), 2:7]
as.data.frame(acomb.true.tax40)

# B_asv78  Actinobacteria     Actinobacteria Streptomycetales Streptomycetaceae  Streptomyces Unclassified
# B_asv91  Actinobacteria     Actinobacteria Streptomycetales Streptomycetaceae  Streptomyces Unclassified
# B_asv93  Actinobacteria     Actinobacteria Streptomycetales Streptomycetaceae  Streptomyces Unclassified
# B_asv110 Actinobacteria     Actinobacteria Streptomycetales Streptomycetaceae  Streptomyces Unclassified
# B_asv126 Actinobacteria     Actinobacteria Streptomycetales Streptomycetaceae  Streptomyces Unclassified
# B_asv128 Actinobacteria     Actinobacteria Streptomycetales Streptomycetaceae  Streptomyces Unclassified
# B_asv165 Proteobacteria Betaproteobacteria  Methylophilales  Methylophilaceae Methylotenera Unclassified

#### Calculate which samples have more/less of the different taxa
view(out.40$feature_table)
feat.flip <- t(out.40$feature_table) ## flip feature table
meta.filt <-nz.16s.Phylo.age@sam_data[which(rownames(nz.16s.Phylo.age@sam_data) %in% rownames(feat.flip)),]
meta.filt 

taxa.tab <- cbind(feat.flip, meta.filt) ## merge metadata and feature table

young.tab.tax <- subset(taxa.tab, YoungvsOld %in% c("Young"))
old.tab.tax <- subset(taxa.tab, YoungvsOld %in% c("Old"))

mean(young.tab.tax$B_asv78) ##13.1875
mean(old.tab.tax$B_asv78) ## 231.6667

mean(young.tab.tax$B_asv91) ## 15.9375
mean(old.tab.tax$B_asv91) ## 211.25

mean(young.tab.tax$B_asv93) ## 18.875
mean(old.tab.tax$B_asv93)## 189.6667

mean(young.tab.tax$B_asv110) ## 12.25
mean(old.tab.tax$B_asv110) ## 163.4167

mean(young.tab.tax$B_asv126) ## 10.75
mean(old.tab.tax$B_asv126)## 116

mean(young.tab.tax$B_asv128) ## 10.75
mean(old.tab.tax$B_asv128) ## 102.5833

mean(young.tab.tax$B_asv165) ## 42.75
mean(old.tab.tax$B_asv165)## 3.25

################## CREATE A HEATMAP OF ASV DATA
###### REORGANIZE DATA FOR HEATMAP
## need to average the samples by young vs. old 
taxa.tab2 <- taxa.tab[,-c(66:76)]
avgs <- aggregate(. ~ YoungvsOld, data=taxa.tab2, FUN= mean)
str(avgs)
row.names(avgs) <- avgs$YoungvsOld
avgs <- avgs[,-1]

## make dataframe numeric
avgs2 <- avgs %>% mutate_at(1:20, as.numeric)
str(avgs2) ## sums2 now numeric
avgs2 <- as.matrix(avgs2)

avgs2 <- t(avgs2)

## change ASVs to genus/species
heat.tax <-red.taxon[which(rownames(red.taxon) %in% rownames(avgs2)), 6:7]
heat.tax$con <- paste(heat.tax$Genus, heat.tax$Species)

# ##rename the unclassified, asvs with same species label, and indicate * for sig diff
heat.tax[34, 3] = "Streptomyces Unclassified *"
heat.tax[45, 3] = "Streptomyces Unclassified *"
heat.tax[46, 3] = "Streptomyces Unclassified *"
heat.tax[58, 3] = "Streptomyces Unclassified *"
heat.tax[62, 3] = "Streptomyces Unclassified *"
heat.tax[63, 3] = "Streptomyces Unclassified *"
heat.tax[65, 3] = "Methylotenera Unclassified *"


row.names(avgs2) <- heat.tax$con
view(avgs2)
rownames(avgs2)

col_order <- c("Young", "Old")
avgs2 <- avgs2[, col_order] ## re order young vs. old

## order rows alphabetical
rr <- avgs2[order(row.names(avgs2)), ]

## heatmap (> 25% abundance)
pheatmap(rr, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = hcl.colors(100, "Tofino"))

###################### ANCOMBC: Phylum Level
# Aggregate to phylum level
phylum_data = aggregate_taxa(nz.16s.Phylo.age, "Phylum")
# The taxonomy table
tax_mat = as(tax_table(phylum_data), "matrix")

### Step 2: run ancombc function: zero_cut threshold is "0.99" 
### so includes all samples/taxa

out.phy = ancombc(phyloseq = phylum_data, formula = "YoungvsOld",
                  p_adj_method = "holm", zero_cut = 0.99, lib_cut = 0,
                  group = "YoungvsOld", struc_zero = TRUE, neg_lb = FALSE,
                  tol = 1e-5, max_iter = 100, conserve = TRUE,
                  alpha = 0.05, global = FALSE)
resphy = out.phy$res
tab_diffphy = resphy$diff_abn
tab_phy = resphy$p_val
tab_phy = resphy$q_val ## some significance

as.data.frame(tab_diffphy)
count(tab_diffphy, tab_diffphy$YoungvsOld == "TRUE")
true.classphy <- subset(tab_diffphy, tab_diffphy$YoungvsOldYoung == "TRUE")
as.data.frame(true.classphy)


### Step 2: run ancombc function again, with zero_cut threshold "0.75"
### so only includes taxa in > 25% of samples

out.phy2 = ancombc(phyloseq = phylum_data, formula = "YoungvsOld",
                   p_adj_method = "holm", zero_cut = 0.75, lib_cut = 0,
                   group = "YoungvsOld", struc_zero = TRUE, neg_lb = FALSE,
                   tol = 1e-5, max_iter = 100, conserve = TRUE,
                   alpha = 0.05, global = FALSE)
resphy2 = out.phy2$res
tab_diffphy2 = resphy2$diff_abn
tab_phy2 = resphy2$p_val
tab_phy2 = resphy2$q_val ## some significance

as.data.frame(tab_diffphy2)
count(tab_diffphy2, tab_diffphy2$YoungvsOld == "TRUE")
true.classphy2 <- subset(tab_diffphy2, tab_diffphy2$YoungvsOldYoung == "TRUE")
as.data.frame(true.classphy2)

