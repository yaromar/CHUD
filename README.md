# CHUD

AnnotatingVariants.py
Developer: Maryam Zekavat
Function: Used this to annotate a list of 12.5 million unique somatic variants in Hail-0.2. Script is meant to be run on the google cloud.
Input: List of unique variants
Output: List of unique variants with annotations

Merge_SomaticVCFS_1_10000.py
Developer: Maryam Zekavat
Function: Used this to convert somatic variants from within 50,000 variant call format (VCF) files to tab-delimited .bgz files for further analysis. Additionally, this script also performs quality control filtration of the somatic variants to high-quality variants.
Input: 50,000 vcf.gz files
Output: 50,000 .bgz files (tab-delimited)

Rscript_regression.R
Developer: Maryam Zekavat
Function: Used this to perform associations of variant counts among various groupings with age, CHIP clone size, and phenotypes (AML, MPN, CAD). This script was also used for the plots created in the paper.
Input: Phenotype file and variant counts
Output: Logistic regression results, CoxPH association results with CAD, and figures

Chuddb_init_script.sql
Developer: Vimig 
Function:
Input:
Output:

Data_exploration.py
Developer: Yaro
Function:
Input:
Output:

variantCount.py
Developer: Yaro
Function:
Input:
Output:

variantCountBatch.py
Developer: Vimig
Function:
Input:
Output:

Vis_clustering_yaro.ipynb
Developer: Yaro
Function:
Input:
Output:
