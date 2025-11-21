### Differential Expression in R

## Libraries
library(DESeq2)
library(tidyverse)

## Importing files for DE analysis
# salmon read counts
txi <- readRDS("RObjects/txi.rds")
head(txi$counts)
# sample sheet
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", col_types = 'cccc')

all(colnames(txi$counts) == sampleinfo$SampleName)

## Simple Model 
simple.model <- as.formula(~TimePoint)
model.matrix(simple.model, data=sampleinfo)

simple.model <- as.formula(~Status)
model.matrix(simple.model, data=sampleinfo)

# relevel 'Infection' in sampleinfo
sampleinfo <- mutate(sampleinfo, Status=fct_relevel(Status, "Uninfected"))
model.matrix(simple.model, data=sampleinfo)


## Construct DESeq Dataset
### design model, raw counts, metadata

ddsObj.raw <- DESeqDataSetFromTximport(txi=txi, colData=sampleinfo, design = simple.model)

### filtering
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep, ]

## Perform Diff Exp Analysis
# DESeq() 

# estimate size factors
ddsObj <- estimateSizeFactors(ddsObj.filt)
normalizationFactors(ddsObj)

# visualize the effect of pre/post 'normalization'
logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)
limma::plotMA(logcounts, array = 5, ylim=c(-5,5))
abline(h=0, col='red')

logcounts <- log2(counts(ddsObj, normalized = TRUE) + 1)
limma::plotMA(logcounts, array = 5, ylim=c(-5,5))
abline(h=0, col='red')

# estimate gene dispersion
ddsObj <- estimateDispersions(ddsObj)
plotDispEsts(ddsObj)

# fit negative binomial distribution
ddsObj <- nbinomWaldTest(ddsObj)

### Results
results.simple <- results(ddsObj, alpha=0.05)
results.simple
# upregulated
sum(results.simple$log2FoldChange>0 & results.simple$padj < 0.05, na.rm = TRUE)
# downregulated
sum(results.simple$log2FoldChange<0 & results.simple$padj < 0.05, na.rm = TRUE)


## Additive model (2 factors)
additive.model <- as.formula(~TimePoint + Status)

ddsObj.raw <- DESeqDataSetFromTximport(txi=txi, colData=sampleinfo, design=additive.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep, ]

### Exercise 3
ddsObj <- DESeq(ddsObj.filt)
res.additive <- results(ddsObj, alpha = 0.05)

# how many coefficients? 2
model.matrix(additive.model, data=sampleinfo)
# what is the reference group? d11, uninfected

# what contrasts?
# Timepoint (d11 vs d33)
# Infection (uninfected vs infected)
res.additive
# num padj < 0.05
sum(res.additive$padj < 0.05, na.rm = TRUE)

# look at different contrasts
resultsNames(ddsObj)

results.InfectedvUninfected <- res.additive
rm(res.additive)

# order genes on padj
topGenesIvU <- as.data.frame(results.InfectedvUninfected) %>%
  rownames_to_column("GeneID") %>% 
  top_n(100, wt = -padj)

# Exercise 5
resultsNames(ddsObj)

# pull out results for time point

res_d33_d11 <- results(ddsObj, alpha=0.05, name = "TimePoint_d33_vs_d11")
res_d33_d11

sum(res_d33_d11$padj < 0.05, na.rm=TRUE)

### PCA /rlog
vstcounts <- vst(ddsObj.raw, blind=TRUE)
plotPCA(vstcounts, intgroup = c("Status", "TimePoint"))

### Exercise 6
int.model <- as.formula(~TimePoint * Status)
ddsObj.raw <- DESeqDataSetFromTximport(txi=txi, colData=sampleinfo, design = int.model)
keep <- rowSums(counts(ddsObj.raw)) > 5 
ddsObj.filt <- ddsObj.raw[keep, ]

ddsObj.interaction <- DESeq(ddsObj.filt)

results.int <- results(ddsObj.interaction, alpha=0.05)
results.int

resultsNames(ddsObj.interaction)

# what genes are expressed between uninfected and infected mice at Day 11
results.interaction.11 <- results(ddsObj.interaction, 
                                  alpha=0.05,
                                  name = "Status_Infected_vs_Uninfected")

# what genes are expressed between uninfected and infected mice at Day 33
results.interaction.33 <- results(ddsObj.interaction, 
                                  alpha=0.05,
                                  contrast = list(c("Status_Infected_vs_Uninfected", "TimePointd33.StatusInfected")))


sum(results.interaction.11$padj < 0.05, na.rm=TRUE)
sum(results.interaction.33$padj < 0.05, na.rm=TRUE)

### Exercise 7
resultsNames(ddsObj.interaction)

# day 11 vs d33, infected mice
results.infected.d33d11 <- results(ddsObj.interaction,
                                   alpha=0.05,
                                   contrast = list(c("TimePoint_d33_vs_d11", "TimePointd33.StatusInfected")))

# day 11 vs d33, uninfected mice
results.uninfected.d33d11 <- results(ddsObj.interaction,
                                   alpha=0.05,
                                   name = "TimePoint_d33_vs_d11")

sum(results.uninfected.d33d11$padj < 0.05, na.rm=TRUE)

sum(results.infected.d33d11$padj < 0.05, na.rm=TRUE)

# save results
saveRDS(ddsObj.interaction, "DESeqDataSet.interaction.rds")
saveRDS(results.interaction.11, "DESeqResults.interaction_d11.rds")
saveRDS(results.interaction.33, "DESeqResults.interaction_d33.rds")


