# CHUD Project

Welcome to the working directory of the CHUD Project! The goal of this project is to identify new somatic mutations that may influence clonal hematopoiesis of indeterminate potential (CHIP). These novel mutations are sometimes called CHUD mutations. What follows is a description of the files by the authors and links to our papers for the Yale CBB 750 Course: Core Topics in Biomedical Informatics (Spring 2020).

## Links to accompanying documents

**Final Paper (must be Yale affiliate):**  https://docs.google.com/document/d/1iSKAuRb69XYnhAc-ekac-dC_p3VXT_xXHuklum3oSbs/edit?usp=sharing   
**Final Presentation (must be Yale affiliate):** https://docs.google.com/presentation/d/1RdI-UUncxgI_9gh3H0MrojjGzQpWJ6WO9-MKAtkrKBo/edit?usp=sharing    
**Link to Github project:** https://github.com/yaromar/CHUD    

------

## Installing the code

We tried to make it as easy as possible!  

1. Make sure you have [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installed 
2. Create a new environment using the command `conda env create -f chudenv.yml` from within your git directory.
3. Activate your new conda env using `conda activate myenv` 
4. Verify it worked using `conda env list`    
5. Done!   

-----

## Description of files

### Data Manipulation Scripts 

_AnnotatingVariants.py_
Developer: Maryam Zekavat  
Function: Used this to annotate a list of 12.5 million unique somatic variants in Hail-0.2. Script is meant to be run on the google cloud.  
Input: List of unique variants  
Output: List of unique variants with annotations    

_Merge_SomaticVCFS_1_10000.py_  
Developer: Maryam Zekavat  
Function: Used this to convert somatic variants from within 50,000 variant call format (VCF) files to tab-delimited .bgz files for further analysis. Additionally, this script also performs quality control filtration of the somatic variants to high-quality variants.    
Input: 50,000 vcf.gz files  
Output: 50,000 .bgz files (tab-delimited)   

_Chuddb_init_script.sql_  
Developer: Vimig   
Function:   
Input:  
Output:    

_Socrates_CHUD_Data_Cleaner.ipynb_     
Developer: Vimig   
Function:   
Input:  
Output:    

### Variant Count and Analysis Scripts 

_variantCount.py_   
Developer: Yaro  
Function: Filters variants by annotations and produces counts that are used for our regression analyses.  
Input:   
annotation file= 
‘UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz’
        variant calls folder = 
            Filtered_SomaticCalls_v1/
Output:  
variantCount_output.txt    

_variantCountBatch.py_   
Developer: Vimig  
Function:  
Input:  
Output:   

_Rscript_regression.R_  
Developer: Maryam Zekavat  
Function: Used this to perform associations of variant counts among various groupings with age, CHIP clone size, and phenotypes (AML, MPN, CAD). This script was also used for the plots created in the paper (Fig 1-6).
Input: Phenotype file and variant counts  
Output: Logistic regression results, CoxPH association results with CAD, and figures   

## Visualization, Dimensionality Reduction, Clustering

_Vis_clustering_yaro.ipynb_  
Developer: Yaro  
Function: data processing, PCA, embedding(t-SNE, UMAP, t-SNE), clustering(K-means, agglomerative, DBSCAN)
Input: training data=’sample_var_phenos_leuk_topmed_v2.feather’  
Output: figures/stats  

## Classification Methods

_CHUD_Data_Model_Trainer.ipynb_     
Developer: Vimig   
Function:   
Input:  
Output:    

_CHUD_multi_model_run.py_    
Developer: Vimig   
Function:   
Input:  
Output:    


_Data_exploration.py_  
Developer: Vimig  
Function:  
Input:  
Output:    

## Other Files

_socrates_chud_jupyter_job.slurm_ - slurm script to start Vimig's Jupyter notebook on the HPC   
_variant_count_batch_job.slurm_ - different way to run job array on HPC     
_variantcountbatchjob.txt_ - jobfile to run job array on HPC     
_yaro_chud_jupyter_job.slurm_ - slurm script to start Yaro's Jupyter notebook on the HPC


