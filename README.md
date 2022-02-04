# ESCA

# ESCAilluminaPKG
ESCA (Easy-to-use SARS-CoV-2 genome Assembler) is built to facilitate the analysis of NGS sequencing data obtained from amplicon-based sequencing methods.
Most dependencies that are necessary to run ESCA package are packages widely used in bioinformatics (Biopython, Samtools, bwa, etc.).
However, we suggest to use the conda environment (available in the "INPUTDATA" directory) to setup all depedencies needed to run ESCA. To 

# INSTALLATION:
0) download and unzip ESCA-main.zip file, wherever you prefer in your PC
1) If you have conda environment management system already installed in your system, open a new shell into INPUTDATA directory and digit
      $ conda env create -f condaESCAenvironment.yml
      $ conda activate escaenv
   Otherwise, if you don't have conda already installed in you PC let's see at :
   https://conda.io/projects/conda/en/latest/user-guide/install/index.html
   to learn more about conda installation and usage.

# ANALYSIS:
1) put ONLY Illumina paired-end files in fastq.gz in the "INPUTDATA" directory
3) open a shell and go the path of "ESCA-main"
4) run the command:
      $bash ./ESCArun

5) find the assembled SARS-CoV-2 genomes, aligned with Wuhan-Wu-1 reference genome,
in the most recent "REULTS" folder, with file name SARScov2consensusALLPlusRefaln.fa

# HINTS: 
1) remember to delete all input files after the end of the analysis
2) control if ESCA environment is active every time you run ESCA
3) in case of problems, please contact cesare.gruber@inmi.it 

