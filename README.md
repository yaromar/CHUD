# CHUD Project

Welcome to the working directory of the CHUD Project! The goal of this project is to identify new somatic variants that may influence clonal hematopoiesis of indeterminate potential (CHIP). These novel variants may be called clonal hematopoiesis of uknown drivers (CHUD). What follows is a description of the files by the authors and links to our papers for the Yale CBB 750 Course: Core Topics in Biomedical Informatics (Spring 2020).

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

**AnnotatingVariants.py**       
_Developer:_ Maryam Zekavat        
_Function:_ Used this to annotate a list of 12.5 million unique somatic variants in Hail-0.2. Script is meant to be run on the google cloud.     
_Input:_ List of unique variants     
_Output:_ List of unique variants with annotations       

**QC_VCFS_1_10000.py**    
_Developer:_ Maryam Zekavat      
_Function:_ Used this to convert somatic variants from within 50,000 variant call format (VCF) files to tab-delimited .bgz files for further analysis. Additionally, this script also performs quality control filtration of the somatic variants to high-quality variants.       
_Input:_ 50,000 vcf.gz files    
_Output:_ 50,000 .bgz files (tab-delimited)     

**Chuddb_init_script.sql**    
_Developer:_ Vimig      
_Function:_ Create the database schema in PostgreSQL    
_Input:_ None   
Output: New PostgreSQL database when run with psql    

**Socrates_CHUD_Data_Cleaner.ipynb**        
_Developer:_ Vimig      
_Function:_ Merge sample variant and phenotype information to create ML data matrix.       
_Input:_  Sample file, annotation file, and phenotype files.     
_Output:_ three binary [feather](https://github.com/wesm/feather) files, for all variants, somatic, leukemia, and TOPMED variants, and just leukemia and TOPMED variants.       

### Variant Count and Analysis Scripts    

**variantCount.py**      
_Developer_: Yaro     
_Function_: Filters variants by annotations and produces counts that are used for our regression analyses.     
_Input_:    
- annotation file = 'UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz'   
- variant calls folder = Filtered_SomaticCalls_v1/      

_Output_: variantCount_output.txt         

**variantCountBatch.py**      
_Developer:_ Vimig     
_Function:_ does same thing as `variantCount.py` but in specified batched     
_Input:_  intput directory, output directory, and how many batches as command line arguments   
_Output:_ variantCount_output[1-numBatches].txt    

**Rscript_regression.R**     
_Developer:_ Maryam Zekavat     
_Function:_ Used this to perform associations of variant counts among various groupings with age, CHIP clone size, and phenotypes (AML, MPN, CAD). This script was also used for the plots created in the paper (Fig 1-6).     
_Input:_ Phenotype file and variant counts        
_Output:_ Logistic regression results, CoxPH association results with CAD, and figures      

### Visualization, Dimensionality Reduction, Clustering    

**Vis_clustering_yaro.ipynb**     
_Developer:_ Yaro     
_Function:_ data processing, PCA, embedding(t-SNE, UMAP, t-SNE), clustering(K-means, agglomerative, DBSCAN)     
_Input:_ training data=’sample_var_phenos_leuk_topmed_v2.feather’     
_Output:_ figures/stats     

### Classification Methods   

**CHUD_Data_Model_Trainer.ipynb**        
_Developer:_ Vimig     
_Function:_ Trains all four models on specified phenotypes, after some initial data manipulation (creating dummy variables, dropping columns) with lots of output and prints some graphs   
_Input:_ training data=’sample_var_phenos_leuk_topmed_v2.feather’     
_Output:_ figures/stats     

**CHUD_multi_model_run.py**       
_Developer:_ Vimig     
_Function:_ Does the same thing as `_CHUD_Data_Model_Trainer.ipynb_` but way cleaner, as a loop over all specified phenotypes.  
_Input:_ training data=’sample_var_phenos_leuk_topmed_v2.feather’     
_Output:_ pretty ROC curves and confusion matrices, along with saved copies of trained models and CSVs of    importance/coefficients. Also saves summary file of accuracy, precision, recall, and F1 scores for all models trained    

### HPC Job Files   

_socrates_chud_jupyter_job.slurm_ - slurm script to start Vimig's Jupyter notebook on the HPC      
_variant_count_batch_job.slurm_ - different way to run job array on HPC         
_variantcountbatchjob.txt_ - jobfile to run job array on HPC        
_yaro_chud_jupyter_job.slurm_ - slurm script to start Yaro's Jupyter notebook on the HPC   


