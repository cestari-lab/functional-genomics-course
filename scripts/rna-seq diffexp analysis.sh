# RNA-seq analysis: differential gene expression analysis
# We will start with sorted.bam files from RNA-seq performed by Nanopore sequencing.
# The steps covered here include read counting using featureCounts and differential expression analysis using edgeR.

# The dataset includes 3 biological replicates of control and treated samples. Ref. Touray et al. 2023, eLife (https://doi.org/10.7554/eLife.89331.3). 

# First, let's install the package subread. 

conda install -n bio_tools -c bioconda subread

# Create output directory
mkdir -p results/rna_seq

# Activate conda environment
conda activate bio_tools

# Read counting using featureCounts:
featureCounts -h # to check help and options

#featureCounts counts mapped reads for genomic features such as genes, exons, etc. 
#The script below will generate the number of read counts for each exon. It will analyze all sorted.bam files.

featureCounts -LO -a genome/Tb427annotation.gtf -o results/rna_seq/counts.txt -F GTF -t exon -g gene_id --readExtension3 200 --fraction  -T4 -Q1 $(ls dataset/rnaseq/*sorted.bam)

# The script below will combine all counts data into a matrix file for edgeR analysis
grep -v "#" results/rna_seq/counts.txt | cut -d$'\t' -f1,7- > results/rna_seq/counts.matrix

#-----------------------------------------------------------

# RNA-seq analysis using edgeR:
# Open a terminal for R, or open RStudio. 

# You may need to install the package edgeR (Limma should be installed with edgeR) first if not installed yet:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

#Load required libraries
library (limma)
library (edgeR)

# Open matrix file generated from featureCounts:
data <- read.table("results/rna_seq/counts.matrix", header = TRUE, sep = "", skip = 0, row.names = 1)

# check the data
head(data)

#change headings of data (make sure the label matches your samples):
names(data)<- c("Mut-Bio1", "Mut-Bio2", "Mut-Bio3", "Wt-Bio1", "Wt-Bio2", "Wt-Bio3")

#Create DGEList and statistical design
y <- DGEList(counts=data, group=group)

# Check number of rows in the count matrix
nrow (y$counts)

# Define the experimental design
group <- factor(c(2,2,2,1,1,1)) # 1=WT, 2=Mutant

# Filter out lowly expressed genes using the following commands:
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Check number of rows in the count matrix - note any difference in the number of rows?
nrow (y$counts)

# Define the experimental design
design <- model.matrix(~group) # 1=WT, 2=Mutant, comparison will be Mutant vs WT (2 vs 1)

# Normalization of data
y <- normLibSizes(y)

#Generate plot MSD (Multi-Dimensional Scaling) for sample comparisons
plotMDS.pdf <-plotMDS(y, col=c(rep("black",3), rep("red",3)))

# Save the MDS plot
png("results/rna_seq/MDS_plot.png") 
plotMDS.pdf <-plotMDS(y, col=c(rep("black",3), rep("red",3)))
dev.off()   

# Estimate sample dispersion
y <- estimateDisp(y, design)

#Plot Biological correlate of variation
plotBCV(y)

# Save the BCV plot 
png("results/rna_seq/BCV_plot.png")
plotBCV(y)
dev.off()

# Generate statistical analysis using generalized linear model
fit <- glmQLFit(y, design)

# Compare groups treated vs non-treated
qlf.2vs1 <- glmQLFTest(fit, coef=2) 

# Check the table of differentially expressed genes
head (qlf.2vs1$table)

# Check top differentially expressed genes
topTags(qlf.2vs1)

# write the results to a file
all.results <- as.data.frame(qlf.2vs1$table)    
write.table(all.results, file="results/rna_seq/dge_results_all.txt", sep="\t", quote=FALSE, col.names=NA)

# Plot volcano plot for DEGs
# Install EnhancedVolcano package if not installed yet
BiocManager::install('EnhancedVolcano')

# Load the library
library(EnhancedVolcano)    

EnhancedVolcano(all.results, lab = rownames(all.results), 
    x = "logFC", y = "PValue", title = "Volcano Plot of Differentially Expressed Genes", 
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 3.0,
    labSize = 6.0)

# Save the volcano plot
png("results/rna_seq/volcano_plot.png")
EnhancedVolcano(all.results, lab = rownames(all.results), 
    x = "logFC", y = "PValue", title = "Volcano Plot of Differentially Expressed Genes", 
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 1.0,
    labSize = 3.0)
dev.off()       

#-----------------------------------------------------------
# End of RNA-seq analysis
#-----------------------------------------------------------

