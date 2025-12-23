# Denovo assembly of E. coli genome from Illumina single-end reads.
# We will use the reads we trimmed and filtered in the previous class.
# The assembly will be performed using SPAdes. 
# The quality of the assembly will be assessed using QUAST and BUSCO.
# We will do step-by-step instructions in this script for a better understanding of the process.    

# First let's install SPAdes. We will use conda, as done for the other tools.
# Make sure you have conda installed and the bio_tools environment created as shown in previous scripts.
# You can read more about SPAdes here: https://github.com/ablab/spades
conda install -n bio_tools -c bioconda spades

# Run a test to make sure SPAdes is installed correctly:
spades.py --test

# Check the output files to see if the test ran successfully. You should have a folder called spades_test in your current directory.

# Now we can proceed with the assembly.
# Set the path to your trimmed and filtered reads. This will simplify the command. Be sure to adjust the path accordingly. 
READS="mypath/dataset/ecolireads-trimmed.fastq"

# Set the output directory for the assembly. Do not forget to change 'mypath' to your actual path.
OUTPUT_DIR="mypath/results"

# Run SPAdes help to check the paprameters 
spades.py --help 

# Run SPAdes for single-end reads (Using a 7M reads file this took 50 minutes with 6 threads and 16GB of RAM) 
# Running with 2M reads and 4 threads and 8GB of RAM took around 20 minutes. 
spades.py --s1 $READS -o $OUTPUT_DIR/spades_assembly --careful --threads 4 --memory 8

# We can check the assembly output files in the specified output directory.
# The key output files are the scaffolds.fasta and contigs.fasta.
# The main assembly file is contigs.fasta located in the output directory. You can view the contigs using:
less $OUTPUT_DIR/spades_assembly/contigs.fasta

# We can also check the assembly statistics using QUAST and BUSCO as shown in the previous script.
# Let's start using QUAST:
mkdir -p $OUTPUT_DIR/quast_output 
quast $OUTPUT_DIR/spades_assembly/contigs.fasta -o $OUTPUT_DIR/quast_output --threads 4

# Use reference-based evaluation with QUAST (Optional)
# You can also use QUAST to perform a reference-based evaluation of your assembly.
# Make sure you have the reference genome file (Ecoli_atcc_genome.fna)
# Run QUAST with reference genome   
# quast $OUTPUT_DIR/contigs.fasta -r genome/Ecoli_atcc_genome.fna -o quast_ref_output --threads 4 

# Now let's run BUSCO to assess genome completeness:
# Let's create an output directory for BUSCO results
mkdir -p $OUTPUT_DIR/busco_output
# Because we are working with a bacterial genome. We will use the bacteria_odb10 lineage dataset.
busco -i $OUTPUT_DIR/spades_assembly/contigs.fasta \
      -o results/busco_output \
      -m genome \
      -l bacteria_odb10 \
      -c 4 \
      -f 

# After the runs are complete, you can check the QUAST and BUSCO output files for assembly quality metrics and completeness scores.
# QUAST results will be in the quast_output directory, and BUSCO results will be in the busco_output directory.
# You can open the summary files to see the results:
less $OUTPUT_DIR/quast_output/report.txt
less $OUTPUT_DIR/busco_output/short_summary.specific.bacteria_odb10.busco_output.txt    

# This concludes the genome assembly tutorial using SPAdes, along with quality assessment using QUAST and BUSCO.

# -----------------------------------------------
# Read mapping and coverage

# First, we need to index the assembled genome.
# We will use bwa for mapping reads to the assembled genome.
# Let's create an output directory for the mapping results
mkdir -p $OUTPUT_DIR/mapping_output
# Index the contigs.fasta file
bwa index $OUTPUT_DIR/spades_assembly/contigs.fasta
# Map the trimmed reads to the indexed genome. 
# Important: Do not forget to change the reads to be mapped. You can use your group trimmed reads file. 
bwa mem $OUTPUT_DIR/spades_assembly/contigs.fasta $READS > $OUTPUT_DIR/mapping_output/mapped_reads.sam

# Convert SAM to BAM format
samtools view -bS $OUTPUT_DIR/mapping_output/mapped_reads.sam > $OUTPUT_DIR/mapping_output/mapped_reads.bam
# Sort the BAM file
samtools sort $OUTPUT_DIR/mapping_output/mapped_reads.bam -o $OUTPUT_DIR/mapping_output/sorted_mapped_reads.bam
# Index the sorted BAM file
samtools index $OUTPUT_DIR/mapping_output/sorted_mapped_reads.bam
# Calculate mapping statistics
samtools flagstat $OUTPUT_DIR/mapping_output/sorted_mapped_reads.bam > $OUTPUT_DIR/mapping_output/mapping_stats.txt
# Calculate coverage statistics using bedtools
# Verify that bedtools is installed
bedtools --version
# Calculate coverage statistics
bedtools genomecov -ibam $OUTPUT_DIR/mapping_output/sorted_mapped_reads.bam > $OUTPUT_DIR/mapping_output/coverage_stats.txt

# -----------------------------------------------

# Synteny analysis using minimap2
# Synteny analysis can be performed using minimap2 to align the assembled genome to a reference genome. 
# The reference genome is already on your folder. It is called Ecoli_atcc_genome.fna
# You have already installed minimap2. You can verify using:
minimap2 --version
# Check for commands available using:
minimap2 --help

# Let's create an output directory for synteny results
mkdir -p $OUTPUT_DIR/synteny_output
# Perform synteny analysis
minimap2 genome/Ecoli_atcc_genome.fna $OUTPUT_DIR/spades_assembly/contigs.fasta > $OUTPUT_DIR/synteny_output/bacteria_synteny.paf
# The output file bacteria_synteny.paf contains the alignment information.

# You can visualize the synteny using visualization tools. We will do this in R:
###Open R terminal
# Install packages if not already installed. You will need to do this only once.
install.packages("ggplot2")
# Load the library ggplot2. You need to do this every time you open R and want to use ggplot2. 
library(ggplot2)
# Check what is your current working directory
getwd()
# Change the working directory to where the PAF file is located
setwd("results/synteny_output")
#List the files in the directory to confirm
list.files()
# Read the PAF file 
paf = read.delim("bacteria_synteny.paf", header=FALSE)
# Check the first few lines of the PAF data
head(paf)
# Subset the first 5 alignments for plotting (optional, to reduce plot size)
paf_sub = paf[1:5, ]
# Subset the last 5 alignments for plotting (optional, to reduce plot size)
paf_sub2 = paf[(nrow(paf)-4):nrow(paf), ]

# Create a synteny plot
mysyntenyplot = ggplot(data=paf_sub, aes(x=V3, xend=V4, y=V8, yend=V9, color=factor(V1))) + geom_segment() + labs(x="Reference coordinate", y="Query coordinate")
# Create a new plot with the last 5 alignments added by editing the code above.

# Display the plot
print(mysyntenyplot)

# Save the plot as a PNG file
ggsave("bacteria_synteny_plot.png", plot=mysyntenyplot)
