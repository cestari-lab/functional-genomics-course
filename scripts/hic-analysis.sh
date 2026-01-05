# Hi-C data analysis script using HiCExplorer

# We will assume the Hi-C data is already produced and will start with .cool files.
# The first steps (not covered) entail alignment of reads to the reference genome and generation of Hi-C matrices.
# We will cover the theory in class, additional scripts can be found in https://github.com/cestari-lab/Hi-C
# These steps are somewhat similar to those covered in other classes. Also, the dataset is quite large which would take to long for us to process here.
# The data used is a subset from Antunes et al. 2025, Nature Communications (https://doi.org/10.1038/s41467-025-66824-3)

# First, let's install cooler, hicexplorer, and pygenometracks. 
# We will create a new conda environment named hic_tools for this purpose, as its dependencies may conflict with other tools.
conda create -n hic_tools -c bioconda -c conda-forge cooler hicexplorer pygenometracks

# Let's first create output directory
mkdir -p results/hic_analysis

# Let's check the data we have using cooler tool
cooler info dataset/hic-data.cool 

# Step 1: Visualize the Hi-C matrix using HiCExplorer 
# Obtain library information
hicInfo -m dataset/hicdata/hic-data.cool --outFileName results/hic_analysis/hic-data-info.txt

# Plot distance decay plot
hicPlotDistVsCounts -m dataset/hicdata/hic-data.cool -o results/hic_analysis/hic-data-distance-decay.png

# Let's also plot the relationship between short vs long-range contacts. You could repeat this with the normatlized and corrected matrices later.
hicPlotSVL -m dataset/hicdata/hic-data.cool --distance 500000 --threads 6 --plotFileName results/hic_analysis/hic-data_plotSvL.png --outFileName results/hic_analysis/hic-data_SvL_pvalues.txt --outFileNameData results/hic_analysis/hic-data_SvL_rawData.txt 

# Change the matrix resolution. For example, from 10Kb, 50Kb and 100Kb.
# The current matrix is at 1Kb resolution.
# After plot matrices and compare their contact maps at different resolutions.
hicMergeMatrixBins -m dataset/hicdata/hic-data.cool -o dataset/hicdata/hic-data-10k.cool -nb 10     

# Try a diagnostic plot to see the quality of the data, and identify biases before normalization and correction.
hicCorrectMatrix diagnostic_plot --matrix dataset/hicdata/hic-data-10k.cool -o results/hic_analysis/hic-data-10k-diagnostic.png

# Normalize the matrix using the coverage normalization method
hicNormalize -m dataset/hicdata/hic-data-10k.cool --normalize smallest  -o dataset/hicdata/hic-data-10k-norm.cool

# Balance the matrix using Knight-Ruiz (KR) matrix method
hicCorrectMatrix correct -m dataset/hicdata/hic-data-10k-norm.cool -o dataset/hicdata/hic-data-10k-nkr.cool --correctionMethod KR

# Now plot matrices and compare their contact maps at different resolutions, e.g. 10K and 100K.
hicPlotMatrix -m dataset/hicdata/hic-data-10k-nkr.cool --clearMaskedBins -t 'hic-data-10kb-res' --log --colorMap OrRd -o results/hic_analysis/hic-data-10k-nkr.png

# You can change colors and other parameters as needed. For example, to plot only a specific region:Check the chromosome names first
cat dataset/chrom.sizes | grep "Chr*"
# Now plot a specific region, e.g., Chr11
hicPlotMatrix -m dataset/hicdata/hic-data-50k-nkr.cool --clearMaskedBins -t 'hic-data-50kb-res' --log --colorMap OrRd -o results/hic_analysis/hic-data-50k-nkr-chr11.png --region Chr11_5A_core_3A_Tb427v9

#-----------------------------------------------------------

# Finding compartments, TADs, and loops. 

# -----------------------------------------------------------
# Analyze chromosome compartments using HiCExplorer
# Regions in a Hi-C matrix can generally be assigned to either the active or the inactive compartment, also called ‘A’ and ‘B’ compartments, respectively. 
# The eigenvector of the correlation matrix is used to derive compartment type and strength for each matrix bin. 
# Generally, regions with positive values are assigned the ‘A’ and regions with negative values the ‘B’ compartment.  

mkdir -p results/hic_analysis/compartments

# Generate observed/expected matrix and plot it
hicTransform -m dataset/hicdata/hic-data-50k-nkr.cool --method obs_exp -o dataset/hicdata/hic-data-50k-nkr-exp.cool

# Identify chromosome compartments using PCA analysis with HiCExplorer
# First, compute the eigenvectors
hicPCA --matrix dataset/hicdata/hic-data-50k-nkr-exp.cool --whichEigenvectors 1 2  -f bedgraph -o results/hic_analysis/compartments/hic-data-50k-nkr-pca1.bedgraph results/hic_analysis/compartments/hic-data-50k-nkr-pca2.bedgraph
# Repeat to generate a text file with eigenvector values in bigwig format (it will be used for plotting later)
hicPCA --matrix dataset/hicdata/hic-data-50k-nkr-exp.cool --whichEigenvectors 1 2  -f bigwig -o results/hic_analysis/compartments/hic-data-50k-nkr-pca1.bw results/hic_analysis/compartments/hic-data-50k-nkr-pca2.bw

# Then use the PCA to generate the global compartmentalization plot
hicCompartmentalization --obsexp_matrices dataset/hicdata/hic-data-50k-nkr-exp.cool --pca results/hic_analysis/compartments/hic-data-50k-nkr-pca1.bedgraph -o results/hic_analysis/compartments/hic-data-50k-nkr-compartment.png 

# Plot the matrix along with the first eigenvector for a specific chromosome, e.g., Chr8 to observe compartments per chromosome
hicPlotMatrix -m dataset/hicdata/hic-data-50k-nkr.cool --clearMaskedBins -t 'hic-data-50kb-res-comp' --log --vMax 1000 --colorMap OrRd -o results/hic_analysis/compartments/hic-data-50k-nkr-exp-comp-chr10.png --bigwig results/hic_analysis/compartments/hic-data-50k-nkr-pca1.bw --region Chr10_5A_core_3A_Tb427v9    

#-----------------------------------------------------------
# Identify loops and TADs using hicExplorer. 

# Create output directory for loops and TADs
mkdir -p results/hic_analysis/loops_tads

# Find loops. Read the documentation for the parameters
hicDetectLoops -m dataset/hicdata/hic-data-10k-nkr.cool -o results/hic_analysis/loops_tads/hic-data-10k-nkr-loops.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 

#Find TADs. Read the documentation for the parameters
hicFindTADs -m dataset/hicdata/hic-data-10k-nkr.cool --outPrefix results/hic_analysis/loops_tads/hic-data-10k-nkr_tad_thres0.005_delta0.01_fdr --minDepth 30000 --maxDepth 100000 --step 10000 --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr -p 6

# Plot the results
hicPlotMatrix -m dataset/hic-data-50k-nkr.cool --clearMaskedBins -t 'hic-data-50kb-res-comp' --log --vMax 1000 --colorMap Blues -o results/hic_analysis/loops_tads/hic-data-50k-nkr-exp-comp-tads-loops-chr10.png --bigwig results/hic_analysis/compartments/hic-data-50k-nkr-pca1.bw --region Chr10_5A_core_3A_Tb427v9 --loops results/hic_analysis/loops_tads/hic-data-10k-nkr.bedgraph --tads results/hic_analysis/loops_tads/hic-data-10k-nkr_tad_thres0.005_delta0.01_fdr_domains.bed


# Use pyGenomeTracks to plot the results. You need to generate track_triangle.ini. See documentation: https://pygenometracks.readthedocs.io/en/latest/
# I generated a sample .ini file for you to use as a template: dataset/pyGenomeTrack.ini
# Edit this file to point to your data files
nano dataset/pyGenomeTrack.ini

# Check the .ini files for details on how to set up the tracks. Here is an example command to plot loops and TADs along with the Hi-C map for a specific region.
pyGenomeTracks --tracks dataset/hicdata/pyGenomeTrack.ini --region Chr9_5A_core_3A_Tb427v9:1-4020287 -o results/hic_analysis/loops_tads/hic-data-all_chr9.png

# To compare, I included a dataset with TADs using from a richer resolution Hi-C dataset 
# The file is located at dataset/hicdata/TAD_domains.bed 
# Replace TAD file path accordingly in the .ini file and re-run the command above to generate a comparison plot.

# End of script
