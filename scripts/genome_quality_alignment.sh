# Genome Quality Comparison and Alignment

# First, we will compare the quality of a few genomes
# We have two genomes assembled with different parameters or different assemblers. We will compare their quality using QUAST. 
#The genomes are from Trypanosoma cruzi Sylvio X10 strain. One was generated in 2018 and the other in 2025 using different sequencing technologies and assembly methods.
quast TcSylvio_2018_genome.fasta TcSylvio_2025_genome.fasta -o tc_quast_comparison_out

# The output folder 'tc_quast_comparison_out' will contain various reports and plots comparing the two assemblies, including N50, number of contigs, total length, GC content, and more.
# You can visualize the results by opening the 'report.html' file in the output folder with a web browser.
# You can also install the extension "vscode-pdf-viewer" in VS Code to view pdf files directly.
# Now let's have a look at the read quality using FastQC and seqkit
# We will first have a look at a few sequences using seqkit
seqkit head -n 5 ecolireads-group1.fastq.gz 

# Then, get some statistics about the reads
seqkit stats ecolireads-group1.fastq.gz > ecolireads-group1_stats.txt

# Now run FastQC to assess the quality of the reads (this is a more complete analysis)
# Let's create an output folder for FastQC results
mkdir fastqc-group1-output
# Now run FastQC
fastqc ecolireads-group1.fastq.gz -o fastqc-group1-output
# FastQC will generate an HTML report in the output folder 'fastqc-group1-output'. 
# Open the HTML file in a web browser to visualize the quality metrics of the reads.   

# Now we will trim the reads using Trimmomatic
trimmomatic SE -threads 4 \
  ecolireads-group1.fastq.gz \
  ecolireads-group1_trimmed.fastq.gz \
  ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
  SLIDINGWINDOW:4:20 \
  MINLEN:50 

#-----------------------------------------------------------

# Genome Alignment and Quality Assessment
# Now we will align sequencing reads to a reference genome and assess the alignment quality.
# There are several alignment tools available. For example, minimap2, BWA, Bowtie, etc. Here we will use BWA for this purpose: https://github.com/lh3/bwa
# Later on, we will explore other aligners as well. 

# Here we use BWA as an example
bwa index Ecoli_atcc_genome.fna
bwa mem Ecoli_atcc_genome.fna ecolireads-group1_trimmed.fastq.gz > ecoli-group1-aligned_reads.sam

# Then convert SAM to BAM, sort, and index
samtools view -bS ecoli-group1-aligned_reads.sam > ecoli-group1-aligned_reads.bam
samtools sort ecoli-group1-aligned_reads.bam > ecoli-group1-aligned_reads_sorted.bam
samtools index ecoli-group1-aligned_reads_sorted.bam

# Now we can assess the alignment quality using samtools flagstat
samtools flagstat ecoli-group1-aligned_reads_sorted.bam > alignment_stats.txt        

# Now you can filter for secondary and supplementary alignments if needed
# You can check for the flags here: https://www.htslib.org/doc/samtools-flags.html
# Also some useful information here: https://broadinstitute.github.io/picard/explain-flags.html
samtools view -b -F 2308 ecoli-group1-aligned_reads_sorted.bam > ecoli-group1-filtered_aligned.bam
samtools index ecoli-group1-filtered_aligned.bam 

# Repeat the flagstat on the filtered alignments
samtools flagstat ecoli-group1-filtered_aligned.bam > ecoli-group1-filtered_alignment_stats.txt

# Now we can check the 'alignment_stats.txt' file for detailed statistics about the alignment, including total reads, mapped reads, properly paired reads, and more.
# This information will help us evaluate the quality of the genome assembly based on how well the reads align to it.
# You can also calculate coverage statistics,
samtools coverage ecoli-group1-aligned_reads_sorted.bam > ecoli-group1_coverage_stats.txt
# Try it again using the option -m or -A.  

# Review 'coverage_stats.txt' for coverage depth and breadth information.
# This will give you insights into how well the genome is covered by the reads, which is another critical metric for assessing genome quality.

# You can visualize coverage using tools like IGV or generate coverage plots using R or Python based on the coverage statistics obtained.

# We will subset the alignment based on mapping quality if needed
samtools view -b -q 30 ecoli-group1-aligned_reads_sorted.bam > high_quality_aligned.bam
samtools index high_quality_aligned.bam 


# Now you can re-assess the alignment statistics and coverage for the high-quality alignments
samtools flagstat high_quality_aligned.bam > high_quality_alignment_stats.txt
samtools coverage high_quality_aligned.bam > high_quality_coverage_stats.txt    

# Review these files for improved alignment and coverage metrics after filtering for high-quality mappings.
# This concludes the genome quality comparison and alignment assessment.

# You can visualize coverage using tools like IGV or generate coverage plots using R or Python based on the coverage statistics obtained.
# This concludes the genome quality comparison and alignment assessment

