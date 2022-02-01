
# Introduction 

We provide a snakemake pipeline for rRNA modification analysis of Nanopore sequencing reads based on multiple conditions. 

Pairwise condition analysis of rRNA modification is carried out using JACUSA2 tool with call-2 option. The output is a BED file with many information including Mismatch, Insertion, and Deletion scores reflecting variant discrimination. This output file will be used for the downstream analysis of rRNA modification detection and its evaluation.
# Installation
We recommend to install software dependencies via `Conda` on Linux. You can find Miniconda installation instructions for Linux [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
Make sure you install the [Miniconda Python3 distribution](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
For performance and compatibility reasons you should install `Mamba` via conda to install Snakemake. See [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more details.
```
conda install -c conda-forge mamba
```
Once you have installed Conda and Mamba, you can download the Snakemake pipeline and the example datasets.
```
git clone https://github.com/dieterich-lab/rRNA_by_ONT.git
cd rRNA_by_ONT/Computational\ Analysis/
```

Then, you install the required packages after creating an environment with Snakemake installed `environment.yml` file.
```
mamba env create -f environment.yaml
conda activate rRNA_by_ONT_env
```
Before executing the Snakemake workflow, download JACUSA2 [jar](https://github.com/dieterich-lab/JACUSA2) file and make sure that you set the path `path_jar` in the config file.
# General Description

The main functions for the downstream analysis are within notebooks in the `notebooks` folder. Please make sure that all notebooks are within `notebooks` folder.

* `ModelAnalysis.py.ipynb` for the 3-way JACUSA2 call-2 output analysis of IVT, MUT, and WT samples. The output is a tabular result `Features_JACUSA2CALL2.csv` summarizing different combinations of features that could be used to detect the rRNA modifications.


```
Label   Dtype   Ref_Pos                 Ref                 Pos   Coverage  Mis                  Mis+Del+Ins        MisContext+Del+Ins  Mis+Del+Ins_Context ModStatus Status
wt_ivt  MinION  NR_003286_RNA18SN5_542  NR_003286_RNA18SN5  542   104.0     0.03635630767894327  1.337610673650488  3.151885341090008   9.346965268620096   Unm       Unm
wt_ivt  MinION  NR_003286_RNA18SN5_1850 NR_003286_RNA18SN5  1850  201735.0  35.675751981674686   81.78742812431301  198.69201401204919  337.1216663613741   m62A      m62A
wt_ivt  MinION  NR_003286_RNA18SN5_1852 NR_003286_RNA18SN5  1852  306211.0  28.36449200147763    46.573838118813    160.37353663012618  291.0499321196112   Unm       NA 
```

The columns represent:
  1. Label of the pairwise analysis
  2. Label of the sequencing device
  3. Combination of reference of reads and positions
  4. Reference of reads
  5. Position
  6. Min read coverage between conditions
  8. Mismatch score
  9. Sum of Mismatch, Insertion, and Deletion scores
  10. Sum of Mismatch scores of the position and its neighbors (where distance <3) plus the Insertion and Deletion scores of the observed position. 
  11. Sum of Mismatch, Insertion, and Deletion scores within the 5mer context where the observed position is in the center.
  12. ModStatus and Status reflect the modification status of the observed site and its neighbors. It should be provided as a TSV file `mod_status_file` with the following format: 
  
```
Ref                Pos ModStatus Status
NR_003286_RNA18SN5 1   Unm       Unm
NR_003286_RNA18SN5 32  Unm       NA
NR_003286_RNA18SN5 33  Unm       NA
NR_003286_RNA18SN5 34  psU       psU
NR_003286_RNA18SN5 35  Unm       NA

```
  13. Kmer where the observed site is flanked by 2 positions.
  
Note: This file will be used for the annotation of the modified sites and their neighbors. In case the information is not available, please fill the two columns with "Unm" value.
* `BivariateAnalysis.py` for the analysis of individual JACUSA2 call-2 outputs.
* `IVTvsWT.py.ipynb` for the IVT vs WT analysis (Pairwise analysis).
* `QuantitativeAnalysis.py.ipynb` for the analysis of the impact of
  read coverage on the detection of rRNA modifications. The analysis
  is based on LOF scores and the combined features.
* `bib.py` contains functions for the extraction of features in 5mer context from JACUSA call-2 output where the central position is the observed site and to plot results.
  
### Outlier Detection Method
rRNA modification can be detected as outliers in unidimensional/ multidimensional space represented by the combined features. Local Outlier Factor method is used for this purpose. It requires two main parameters to be set: contamination value `LOF_contamenation` and neighborhood size `LOF_neighbors`.
For bivariate analysis, Interquartile Range IQR method could be used also.
# Usage
To run the snakemake pipeline, a set of parameters for each analysis should be defined in the config file. Inputs are the BAM files obtained from the alignement of basecalled reads to a reference transcriptome using Minimap2. The preprocessing of BAM files in the workflow (Generating the MD tag for JACUSA2 call-2) requires a reference transcriptome, so make sure that you set the path to the reference transcriptome `path_ref` in the config file. 

* For the 3-way JACUSA2 call-2 analysis run  `analysis_aggregate` rule. 
```
snakemake --cores all analysis_aggregate
````
* For the 3-way JACUSA2 call-2 analysis with downsampling of BAM files to different levels of read coverage run  `downsampling_aggregate` rule. The coverage level `coverage` and the seed value `seed` should be specified in the config file.
```
snakemake --cores all downsampling_aggregate
````
This will produce 3 directories: `bam` containing preprocessed bam files, `jacusa` containing JACUSA call-2 output, and `analysis` containing all kinds of tabular results and visualizations for each sample. In the case of downsampling, for each read coverage value $Nb specified in the config file `downsampling_param`, a folder is created under a name ending with **.sampled$Nb**, it contains as many folders as the number of seeds `seed` set in the config file, named **DowS$seed**, $seed is the seed value considered, each folder contains all the analysis outputs: tabular results `Features_JACUSA2CALL2.csv`, scatter plots, and LOF scores.

![](./Plots/GeneticModel/MinION/eg.png "h")

* For WT vs IVT analysis:
```
snakemake --cores all IVTvsWT_analysis
````

* For the analysis of different mixtures of WT/MUT samples replacing the WT sample with a defined read coverage run  `mixing_aggregate` rule. The seed value `seed` and the coverage level to be used `sampling_cov` should be set in the config file.
```
snakemake --cores all mixing_aggregate
````
This will produce a folder in the respective sample directory with the prefix **MixS** and values indicating the seed and the fraction of the WT reads. The folder will contain all the analysis outputs.

* For the bivariate analysis with 2 replicates run `bivariate_analysis` rule. Here, the outlier detection method `method` should be specified in the config file. It could be LOF or IQR.
```
snakemake --cores all bivariate_analysis
````
This will produce a folder containing all the analysis outputs: tabular results `Features_JACUSA2CALL2.csv`, bar plots, and LOF scores.

![alt text](./Plots/BiologicalModel/hbiopsies/lof_0.002/bivariate.png "bivariate_analysis")

Edit the provided `config.yaml` file to match your own input files (WT, IVT, MUT), label, LOF parameters, etc.

# Data
All the BAM files of the genetic model used in the paper are provided in `data` folder. Please make sure that all data are within `data` folder.

