# Installing miniconda and setting up conda environments for bioinformatics tools

# Open Visual Studio Code, then terminal, e.g., Ubuntu (WSL), PowerShell, or Command Prompt.
# I recommend VSCode terminal for the best experience. You should know these steps and instructions in myCourses!

# Additional information is here: https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2

# Download and install Miniconda (latest version for Linux)

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

# After installing, close and reopen your terminal application or refresh it by running the following command:
conda init --all    # You need to do this only once. You may need to close and reopen your terminal after this.  

# Initialize conda by running the following command:
# Note you will see (base) in your terminal prompt, indicating that the base conda environment is active.
source ~/miniconda3/bin/activate

# Deactivate conda environment by running:
conda deactivate

#-----------------------------------------------

# whenever you want to use conda, just run:
source ~/miniconda3/bin/activate
# to deactivate, run:
conda deactivate

# Verify the installation by checking the conda version:
conda --version

# You need to add the bioconda and conda-forge channels to your Conda configuration. 
# Many Bioconda packages require conda-forge and should be given the highest priority. 
# You can do this by running the following commands:
conda config --add channels defaults    
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Now you are ready to use Miniconda and install packages as needed!
#It is highly recommended to install tools in a dedicated environment to avoid conflicts with other software or your base system. Create a new environment named bio_tools (you can choose any name) with the desired tools: 
conda create -n bio_tools minimap2 samtools fastqc flye quast bwa bedtools deeptools seqkit trimmomatic busco python=3.8 -y

#Activate the environment:
conda activate bio_tools    

# Check installed tools:
minimap2 --version
samtools --version
fastqc --version
flye --version
quast --version
bwa
bedtools --version  
deeptools --version
seqkit --version
trimmomatic -version
busco --version
python --version

#When you are done working in the environment, you can deactivate it by running:
conda deactivate

# To install new tools in the bio_tools environment, first activate the environment, or run the install command with the -n flag to specify the environment name. 
# For example, to install the tool "cooler" in the bio_tools environment, run:
# conda install -n bio_tools -c conda-forge -c bioconda cooler

# To get help on tools installed via conda, you can use the --help flag with the tool's command. For example:
minimap2 --help

# Remember to always activate the appropriate conda environment before using the installed tools. 
# This ensures that you are using the correct versions of the tools and their dependencies.


