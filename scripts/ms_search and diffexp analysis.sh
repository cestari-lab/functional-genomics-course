# Analysis of proteomics data
# There are many GUI based tools for proteomics data analysis, such as Trans-Proteomics Piepline, MaxQuant, Scaffold, etc.
# We will use command line tools.
# The tools we will use are: Comet (for database search), Mokapot (for post-processing and FDR estimation), thermorawfileparser (for raw file conversion), and MSstats (for differential expression analysis in R). 
# We will also use some home made python scripts for extracting intensities.  
# Make sure you have conda installed on your system. If not, please install Miniconda or Anaconda first. You should have this if you attended the previous classes. 
# You can download and install Miniconda from https://docs.conda.io/en/latest/miniconda.html

#--------------------------------------------------
# We first need to install the necessary tools using conda

# We will create a conda environment for this purpose
conda create -n ms_tools -c bioconda -c conda-forge comet-ms percolator mokapot pyteomics thermorawfileparser -y
# We may not use percolator in this script, but we will install it for your convenience in case you want to try it out later.

# Activate the environment
conda activate ms_tools

# convert raw files to mzML format using ThermoRawFileParser 
# The options -f=2 specifies mzML indexed output format, whereas --excludeExceptionData avoids writing extra exception data that is not needed for most analyses.
thermorawfileparser -i ms-data.raw -b ms-data.mzML --excludeExceptionData

# Run comet command for options
comet

# Prepare comet parameter file
# Edit the parameters according to your experimental design and data characteristics. Save it as comet.params. 
comet -p 

# Then run comet search after editing the parameter file as needed.
# If you saved parameters with a different name, e.g. comet.params.new, then use: comet -Pcomet.params.new ms-data.mzML
# It will generate a pepXML file named ms-data.pep.xml. The .pin file is for percolator, which we are not using in this script. Select to generate a txt file for visualization. 
comet ms-data.mzML  

# Now we can use mokapot for post-processing and FDR estimation. 
# Check mokapot options using:
mokapot --help

# Run mokapot with appropriate options. Make sure to specify the correct decoy prefix used in your database.
# Here is an example command:
mokapot --decoy_prefix DECOY --proteins tbdbase.fasta ms-data.pep.xml 

# Alternatively, you can run percolator if you want to compare results.
# percolator --help
# percolator -X perc-ms-data.xml ms-data.pin 
# you can find some help at https://github.com/percolator/percolator/wiki/Interface

# Now we have to filter the mokapot results based on a desired FDR threshold.
# You can open the mokapot_results.txt file in excel and filter based on q-value as well.
# It can be any desired threshold, e.g., 0.01 for 1% FDR.

# Extract spectral intensities using pyteomics
# We will use a simple Python script to extract intensities from the mzML file. You can also do this using other tools such as TPP.
# First, we will extract the intensities from the mzML file and save them to a CSV file. Then, we will merge those with the mokapot results to get peptide intensities.
# I prepared two python scripts for you: SpectIntensity.py and XIntensity.py
python SpectIntensity.py ms-data.mzML ms-data.csv

# Then, we will map the intensities to the identified peptides using the following script:
python XIntensity.py ms-data.csv mokapot.psms.txt myresults.csv

# Deactivate the conda environment when done
conda deactivate

# Now, with intentisites, we can analyze differential expression using statistical tools such as MSstats in R. 

#--------------------------------------------------
# R script for differential expression analysis using MSstats

# Install MSstats if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MSstats") 

# Load MSstats library
library(MSstats)

# Check you version
packageVersion("MSstats")

# Get help on MSstats functions
?MSstats    
?MSstats::dataProcess
?MSstats::groupComparison

# List files in the directory to confirm the presence of the intensity data
list.files() 

# Read in the peptide intensity data
data <- read.csv("myresults_intensities.csv", header=TRUE, sep=",")

# Check the structure of the data
str(data)
head(data)

# Convert to MSstats format
msstats_data <- MSstats::dataProcess(data)

# Check processed data
head(msstats_data$FeatureLevelData)

# Setup the contrast matrix
comparison <- matrix(c(-1, 1), nrow=1)
colnames(comparison) <- c("control","tumor")
row.names(comparison)<-"tumor-control comparison"

# Perform differential expression analysis
diff_exp_results <- MSstats::groupComparison(contrast.matrix = comparison, data = msstats_data)

# View results
print(diff_exp_results) # You can further visualize and interpret the results as needed.            

##--------------------------------------------------
# MSstats tutorial script (from MSstats website: https://msstats.org)
# Install MSstats if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MSstats") 

# MSstats installation is complete ~~~

# The datasets is from Chang et al. 2012, Mol Cell Proteomics, doi: 10.1074/mcp.M111.014662) 
# The RawData is ready for MSstats (which already add column Condition, BioReplicate, Run, Intensity). 
# You will need to re-format your data to match this format. 

# Run MSstats (start here if you already install MSstats)
# load the library
library(MSstats)

# Here is how to get the help file
??MSstats    
?MSstats::dataProcess
?MSstats::groupComparison

# Check your working directory
getwd()

# Set your working directory to where your data file is located
setwd("your/path")  # Change to your path  

# List files in the directory
list.files()

# Input data
raw<-read.csv("RawData.oc.csv")
# Check the head of the data
head(raw)
head(raw$Intensity)  # Check your condition column
# Check the structure of the data
str(raw)

# Function: dataProcess
# pre-processing data: quality control of MS runs
? dataProcess
quantData<-dataProcess(raw)

# Check the structure of the processed data
str(quantData)

# check processed data and extract feature-level and protein-level data
head(quantData$FeatureLevelData)
head(quantData$FeatureLevelData$GROUP)  
head(quantData$ProteinLevelData)
head(quantData$ProteinLevelData$Protein)

# output QuantData
# Normalize aggregated protein-level data  
write.csv(quantData$ProteinLevelData,file="QuantData_ProteinData.csv")
# Feature-level (peptide-level) data to csv files
write.csv(quantData$FeatureLevelData,file="QuantData_FeatureData.csv")

# Function: dataProcessPlots
# visualization # change to your protein of interest  
? dataProcessPlots
# To print on screen
dataProcessPlots(data=quantData, type="QCPlot", address = FALSE, which.Protein = "CD59") 
dataProcessPlots(data=quantData, type="ProfilePlot", address = FALSE, which.Protein = "CD59")
dataProcessPlots(data=quantData, type="ConditionPlot", address = FALSE, which.Protein = "CD59")

# You may need to adjust the size of the output plot when printing to PDF or other formats
dataProcessPlots(data=quantData, type="QCPlot", which.Protein = "CD59", width=8, height=4)
dataProcessPlots(data=quantData, type="ProfilePlot", which.Protein = "CD59", width=8, height=4)
dataProcessPlots(data=quantData, type="ConditionPlot", which.Protein = "CD59", width=8, height=4)

# Function: groupComparison
# generate testing results of protein inferences across concentrations
?groupComparison
# Check feature-level data
head(quantData$FeatureLevelData)
# Check the levels of your conditions
levels(quantData$FeatureLevelData$GROUP)
groups = levels(quantData$ProteinLevelData$GROUP)
groups
# Define comparison matrix
# Set row names for the comparison matrix, using the levels of your conditions: -1 for control, 1 for treatment
comparison<-matrix(c(-1,1),nrow=1)
comparison
# Get the levels of your conditions from protein-level data
colnames(comparison) <- groups[order(as.numeric(groups))]
comparison
row.names(comparison)<-"tumor-control comparison"
comparison

# Perform group comparison
resultOneComparison<-groupComparison(contrast.matrix=comparison,data=quantData)
# View the head of the comparison result
head (resultOneComparison$ComparisonResult)
# Write comparison result to CSV file
write.csv(resultOneComparison$ComparisonResult, file="OvarianCancer_SignificanceTestingResult.csv")

# Function: groupComparisonPlots
# visualization for testing results
?groupComparisonPlots

# Visualization 1: Volcano plot
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultOneComparison$ComparisonResult,type="VolcanoPlot",address="OC_Ex1_", width=6, height=6)
# (2) Both FDR cutoff = 0.05; fold change cutoff = 1.5
groupComparisonPlots(data=resultOneComparison$ComparisonResult,type="VolcanoPlot",FCcutoff=1.5,ProteinName=FALSE,address="OC_Ex2_", width=6, height=6)
# Visualization 3: Comparison plot
groupComparisonPlots(data=resultOneComparison$ComparisonResult,type="ComparisonPlot",address="OC_Ex1_", width=6, height=6)

# Visualization 4: Heatmap plot not visible for this dataset as it requires more than 2 conditions.
# groupComparisonPlots(data=resultOneComparison$ComparisonResult,type="Heatmap",address="OC_Ex1_", width=6, height=6)
# End of MSstats tutorial script
