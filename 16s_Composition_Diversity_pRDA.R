## Code for: "The age of absorptive roots impacts root-adjacent microbial composition in grapevines"
### This code covers 16s beta composition, alpha diversity, and partial RDA

## Begin
## Data import and cleaning

## Set directory & check files
setwd("/Users/mjp23/Desktop/RootAgePostDada2") ## insert your directory

### Load required packages ###
## install packages as needed

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
library(microbiome)

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
red.taxon

## Skip transpose step

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

rowSums(asv.perc) # 100
colSums(asv.perc)
sum(colSums(asv.perc))

## Split up age and variation roots from metadata and asv.perc files
## one root per age per vine was selected for analysis; need to separate roots classified as "variation"
## computational analyses were performed and reported for "age" roots, not "variation" roots.
asv.fac <- cbind(metadat.red, asv.perc)
asv.age <- subset(asv.fac, Age.testing.or.variation.testing.%in% c("Age")) ## split into just age-related roots
asv.var <- subset(asv.fac, Age.testing.or.variation.testing.%in% c("Variation")) ## split to just variation-related roots
asv.age.perc <- asv.age[,12:18039] ## new ASV
age.met <- asv.age[,1:11] ## new metadata
asv.var.perc <- asv.var[,12:18039] ## new ASV
asv.var.met <- asv.var[,1:11] ## new metadata

## Just age roots data-set, 
## with residual info(e.g.,taxa not in just-age roots) removed from the full data-set
## w/variation roots 
test <- t(asv.age.perc)
rowSums(test)
no_var <- test[rowSums(test)>0,] ## anything with ALL zeroes removed (due to residual info from the "variation" dataset)

no.var.taxa <- red.taxon[which(rownames(red.taxon) %in% rownames(no_var)),] ## taxonomy
asv.novar <- t(no_var) ## asv
meta.novar <- metadat.red[which(rownames(metadat.red) %in% rownames(asv.novar)),] ## metadat

## Set factors as factors before analysis
# For aage roots
meta.novar$Vine <- as.factor(meta.novar$Vine)
meta.novar$YoungvsOld<-as.factor(meta.novar$YoungvsOld)

## create phyloseq object
novar.16s.ASV.age <- otu_table(asv.novar, taxa_are_rows=F)
novar.16s.Meta.age <- sample_data(meta.novar,errorIfNULL=TRUE)
novar.16s.Taxo.age <- tax_table(as.matrix(no.var.taxa, errorIfNULL=TRUE))
novar.16s.Phylo.age <- phyloseq(novar.16s.ASV.age, novar.16s.Meta.age, novar.16s.Taxo.age) 
novar.16s.Phylo.age ## 7297 total taxa in just-age dataset

save.image("/Users/mjp23/Desktop/RootAgePostDada2/novar_16s.Rdata")

saveRDS(novar.16s.Phylo.age, file = "Novar.16s.rds")
phytest <- readRDS("Novar.16s.rds")

## Aggregate to Phylum and Class levels
phylum_data = aggregate_taxa(novar.16s.Phylo.age, "Phylum")
class_data = aggregate_taxa(novar.16s.Phylo.age, "Class")

phylum_data@tax_table
class_data@tax_table

df.t = data.frame(sample_data(novar.16s.Phylo.age))
df.p = data.frame(sample_data(phylum_data))
df.c = data.frame(sample_data(class_data))

## 1. Beta Composition
## a. PERMANOVA
## ASV
bc.tax = distance(novar.16s.Phylo.age, "bray")
set.seed(10)
tax.perm = adonis(bc.tax ~ df.t$YoungvsOld, strata = df.t$Vine, data= df.t, permuations = 999)
tax.perm

## Phylum
bc.phy = distance(phylum_data, "bray")
set.seed(10)
phy.perm = adonis(bc.phy~ YoungvsOld, data= df.p, permutations = 999)
phy.perm

## b. Graphs
## PCoA
##ASV
ASV.ord.pc = ordinate(novar.16s.Phylo.age, "PCoA", "bray")
pc.ASV <- plot_ordination(novar.16s.Phylo.age, ASV.ord.pc, "samples", color="YoungvsOld")
pc.ASV <- pc.ASV + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2)  +
  labs(title = "PCoA 16s: Selection of roots (n = 44)") +
  theme_bw()
print(pc.ASV)

## with vine included
ASV.ord.pc2 = ordinate(novar.16s.Phylo.age, "PCoA", "bray")
pc.ASV2 <- plot_ordination(novar.16s.Phylo.age, ASV.ord.pc, "samples", color="Vine", shape = "YoungvsOld")
pc.ASV2 <- pc.ASV2 + geom_point(size = 5) + 
  labs(title = "PCoA 16s: Selection of roots (n = 44)") +
  theme_bw()
print(pc.ASV2)

## Phylum
phy.ord.pc = ordinate(phylum_data, "PCoA", "bray")
pc.phy <- plot_ordination(phylum_data, phy.ord.pc, "samples", color="YoungvsOld")
pc.phy <- pc.phy + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(pc.phy)

## CAP
## ASV
tax.ord = ordinate(novar.16s.Phylo.age, "CAP", "bray", ~YoungvsOld)
cap.tax <- plot_ordination(novar.16s.Phylo.age, tax.ord, "samples", color="YoungvsOld")
cap.tax <- cap.tax + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(cap.tax)

## ASV- by vine
vine.tax.ord = ordinate(novar.16s.Phylo.age, "CAP", "bray", ~Vine)
vine.cap.tax <- plot_ordination(novar.16s.Phylo.age, vine.tax.ord, "samples", color="Vine")
vine.cap.tax <- vine.cap.tax + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(vine.cap.tax)

## Phylum
phy.ord = ordinate(phylum_data, "CAP", "bray", ~YoungvsOld)
cap.phy <- plot_ordination(phylum_data, phy.ord, "samples", color="YoungvsOld")
cap.phy <- cap.phy + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(cap.phy)

## 2. Beta Dispersion
## Taxonomic
disp.yo.tax <- betadisper(bc.tax,(as.factor(df.t$YoungvsOld)))
disp.yo.tax
anova(disp.yo.tax) ## no diff 0.6499

## Phylum
disp.yo.phy <- betadisper(bc.phy,(as.factor(df.p$YoungvsOld)))
disp.yo.phy
anova(disp.yo.phy) ## no diff 0.8916


### 3. Alpha Diversity
## Alpha diversity measures on rareified but not changed to percentages data
comp.asv.fac <- cbind(metadat.red, asv.r)
comp.asv.age <- subset(comp.asv.fac, Age.testing.or.variation.testing.%in% c("Age")) ## split into just age-related roots
comp.asv.var <- subset(comp.asv.fac, Age.testing.or.variation.testing.%in% c("Variation")) ## split to just variation-related roots
comp.asv.age.noperc <- comp.asv.age[,12:18039] ## new ASV
comp.age.met <- comp.asv.age[,1:11] ## new metadata
comp.asv.var.noperc <- comp.asv.var[,12:18039] ## new ASV
comp.asv.var.met <- comp.asv.var[,1:11] ## new metadata

comp.16s.ASV.age <- otu_table(comp.asv.age.noperc, taxa_are_rows=F)
comp.16s.Meta.age <- sample_data(comp.age.met,errorIfNULL=TRUE)
comp.16s.Taxo.age <- tax_table(as.matrix(red.taxon), errorIfNULL=TRUE)
comp.16s.Phylo.age <- phyloseq(comp.16s.ASV.age, comp.16s.Meta.age,comp.16s.Taxo.age) ## Phyloseq object for just age-selected roots and no percentages
comp.16s.Phylo.age ## tax will be higher because still includes taxa names from variation data-set

## aggregate at phylum and class levels
phylum_data2 = aggregate_taxa(comp.16s.Phylo.age, "Phylum")
class_data2 = aggregate_taxa(comp.16s.Phylo.age, "Class")

phylum_data2@tax_table
class_data2@tax_table

df.t2 = data.frame(sample_data(comp.16s.Phylo.age))
df.p2 = data.frame(sample_data(phylum_data2))
df.c2 = data.frame(sample_data(class_data2))

## Taxonomic level
rich = estimate_richness(comp.16s.Phylo.age, measures=c("Observed", "Shannon"))
meta<-as.data.frame(comp.16s.Phylo.age@sam_data)
richtable<-cbind(meta,rich)

## observed species qq plot 
qqnorm(richtable$Observed, pch = 1, frame = FALSE)
qqline(richtable$Observed, col = "steelblue", lwd = 2)

## shannon diversity qq plot (to check if normally distributed)
qqnorm(richtable$Shannon, pch = 1, frame = FALSE)
qqline(richtable$Shannon, col = "steelblue", lwd = 2)

anova.obs <-aov(rich$Observed~comp.16s.Phylo.age@sam_data$YoungvsOld)
summary(anova.obs)

anova.shan <-aov(rich$Shannon~comp.16s.Phylo.age@sam_data$YoungvsOld)
summary(anova.shan)

young <- subset(richtable, YoungvsOld %in% c("Young"))
old <- subset(richtable, YoungvsOld %in% c("Old"))
mean(young$Observed) ## 244.9565
mean(young$Shannon) ## 4.946596

mean(old$Observed) ## 261.2857
mean(old$Shannon) ## 5.021599

## Phylum level
rich.phy = estimate_richness(phylum_data2, measures=c("Observed", "Shannon"))
meta.phy <-as.data.frame(phylum_data2@sam_data)
richtable.phy <-cbind(meta.phy,rich.phy)

## observed species qq plot
qqnorm(richtable.phy$Observed, pch = 1, frame = FALSE)
qqline(richtable.phy$Observed, col = "steelblue", lwd = 2)

## shannon diversity qq plot (basically normally distributed)
qqnorm(richtable.phy$Shannon, pch = 1, frame = FALSE)
qqline(richtable.phy$Shannon, col = "steelblue", lwd = 2)

anova.obs.phy <-aov(rich.phy$Observed~phylum_data2@sam_data$YoungvsOld) ## phylum-level but not sure I can do that at this level
summary(anova.obs.phy)

anova.shan.phy <-aov(rich.phy$Shannon~phylum_data@sam_data$YoungvsOld)
summary(anova.shan.phy) ## shannon still insignificant

## Class level
rich.class = estimate_richness(class_data2, measures=c("Observed", "Shannon"))
meta.class <-as.data.frame(class_data2@sam_data)
richtable.class <-cbind(meta.class,rich.class)

## observed species qq plot, not great looking - not sure I can do at phylum level
qqnorm(richtable.class$Observed, pch = 1, frame = FALSE)
qqline(richtable.class$Observed, col = "steelblue", lwd = 2)

## shannon diversity qq plot (basically normally distributed)
qqnorm(richtable.class$Shannon, pch = 1, frame = FALSE)
qqline(richtable.class$Shannon, col = "steelblue", lwd = 2)

anova.obs.class <-aov(rich.class$Observed~class_data2@sam_data$YoungvsOld) ## phylum-level but not sure I can do that at this level
summary(anova.obs.class)

anova.shan.class <-aov(rich.class$Shannon~class_data2@sam_data$YoungvsOld)
summary(anova.shan.class) ## shannon still insignificant

### 4. Composition Information
## Young and old dominant taxa
young.tax <- subset_samples(comp.16s.Phylo.age, YoungvsOld %in% c("Young"))
old.tax <- subset_samples(comp.16s.Phylo.age, YoungvsOld %in% c("Old"))

## sort top young
top50.young <- names(sort(taxa_sums(young.tax), decreasing=TRUE))[1:50]
ps.top50.young <- transform_sample_counts(young.tax, function(OTU) OTU/sum(OTU))
ps.top50.young <- prune_taxa(top50.young, ps.top50.young)

## sort bottom young
bottom50.young <- names(sort(taxa_sums(young.tax), decreasing=TRUE))[1:50]
ps.bottom50.young <- transform_sample_counts(young.tax, function(OTU) OTU/sum(OTU))
ps.bottom50.young <- prune_taxa(bottom50.young, ps.bottom50.young)

## what phyla did the top young taxa belong to
topPhyla.young <- tax_table(ps.top50.young)[, "Phylum"]
topPhyla.young.view <- as(topPhyla.young, "vector")
topyoung <- as.data.frame(topPhyla.young.view)

## what phyla did the bottom young taxa belong to
bottomPhyla.young <- tax_table(ps.bottom50.young)[, "Phylum"]
bottomPhyla.young.view <- as(bottomPhyla.young, "vector")
bottomyoung <- as.data.frame(bottomPhyla.young.view)
print(bottomyoung)

## sort top old
top50.old <- names(sort(taxa_sums(old.tax), decreasing=TRUE))[1:50]
ps.top50.old <- transform_sample_counts(young.tax, function(OTU) OTU/sum(OTU))
ps.top50.old <- prune_taxa(top50.old, ps.top50.old)

## what phyla did the top old taxa belong to?
topPhyla.old <- tax_table(ps.top50.old)[, "Phylum"]
topPhyla.old.view <- as(topPhyla.old, "vector")
topold <- as.data.frame(topPhyla.old.view)

### Subset phylum that are relevant for differential abundance
GAl15 <- subset_taxa(comp.16s.Phylo.age, Phylum=="GAL15")
CB <- subset_taxa(comp.16s.Phylo.age, Phylum=="Candidatus_Berkelbacteria")
Fuso <- subset_taxa(comp.16s.Phylo.age, Phylum=="Fusobacteria")
SR1 <- subset_taxa(comp.16s.Phylo.age, Phylum=="SR1_(Absconditabacteria)")
WW3 <- subset_taxa(comp.16s.Phylo.age, Phylum=="WWE3")

diff.tax <- cbind(WW3@otu_table,SR1@otu_table,Fuso@otu_table,CB@otu_table, GAl15@otu_table)
diff.tax <- t(diff.tax)
diff.tax <- as.data.frame(diff.tax)

all.tax <- comp.16s.Phylo.age@otu_table
p.all.tax <- transform_sample_counts(all.tax, function(x) x / sum(x))
p.all.tax <- as.data.frame(p.all.tax)

p.fil.2 <- subset(p.all.tax, rownames(p.all.tax) %in% rownames(diff.tax))


### Partial RDA
counts.ja <- as.data.frame(asv.age.perc) ### Just age roots, percent asvs
meta.rda.ja <-as.data.frame(age.met) ## again just age roots
meta.rda.env.ja <- meta.rda.ja[ , -which(names(meta.rda.ja) %in% c("OLD.ID","Age.testing.or.variation.testing.","AgeCat"))]

as.factor(meta.rda.env.ja$Row)
is.factor(meta.rda.env.ja$Row) 
class(meta.rda.env.ja$Row)

counts.hel.ja <-decostand(counts.ja, "hellinger")
decorana(counts.hel.ja)

ja.rda.p <- rda(counts.hel.ja ~ YoungvsOld + Condition(Vine), meta.rda.env.ja) 
summary(ja.rda.p)
set.seed(1)
anova(ja.rda.p)
(ja.r2.root.p <- RsquareAdj(ja.rda.p)$r.squared) ##  
(ja.r2.root.p.adj <-RsquareAdj(ja.rda.p)$adj.r.squared) ## 
