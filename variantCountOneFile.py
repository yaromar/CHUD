import pandas as pd
from os import listdir
import argparse

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("annotation-file", help="Full path of annotation file")
# Example: 'UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz'

parser.add_argument("output-file", help="Name of output file")
# Example: 'variantCount_output.txt'

parser.add_argument("square", help="display a square of a given number")                                        
args = parser.parse_args()
print("args", args)


#the annotation file should be in the same directory as the script
annotations = pd.read_csv(args["annotation-file"], delimiter='\t', compression='gzip')[['f0', 'SOMATIC', 'LeukemiaGene', 'TOPMed_CHIPVar']].set_index('f0')
annotations = annotations.fillna({'SOMATIC': '0'})

#the FOLDER with the somatic calls files should be in the current directory
#this folder's name should be 'Filtered_SomaticCalls_v1'
#the files should be compressed and have .bgz extension
#there should be no other files
with open(args["output-file"], 'w') as fh:
    fh.write('\t'.join(['sample_id',
                        'num_FilterMutect',
                        'overlap',
                        'overlap_&_LeukemiaGene',
                        'overlap_&_TOPMed_CHIPVar',
                        'overlap_&_SOMATIC',
                        ####
                        'overlap_&_binom',
                        'overlap_&_LeukemiaGene_&_binom',
                        'overlap_&_TOPMed_CHIPVar_&_binom',
                        'overlap_&_SOMATIC_&_binom',
                        ####
                        'overlap_&_VAF',
                        'overlap_&_LeukemiaGene_&_VAF',
                        'overlap_&_TOPMed_CHIPVar_&_VAF',
                        'overlap_&_SOMATIC_&_VAF',
                        ####
                        'overlap_&_binom_&_VAF',
                        'overlap_&_LeukemiaGene_&_binom_&_VAF',
                        'overlap_&_TOPMed_CHIPVar_&_binom_&_VAF',
                        'overlap_&_SOMATIC_&_binom_&_VAF',
                       ]))
    fh.write('\n')

    for file in listdir('Filtered_SomaticCalls_v1'):
        sample = pd.read_csv('Filtered_SomaticCalls_v1/'+file, delimiter='\t', compression='gzip', usecols=["v", "Binomial_Prob", "VAF"]).set_index('v')
        temp = pd.concat([annotations, sample], axis=1, join='inner')

        sampleID, mutectVar = [file.split('.', 4)[i] for i in (0, 3)]
        sampleID = sampleID.split('.', 4)[0].split('-')[0]


        fh.write('\t'.join([sampleID,
                            mutectVar,
                            str(len(pd.concat([annotations, sample], axis=1, join='inner'))), #overlap
                            str(len(temp.loc[temp['LeukemiaGene'] == 1])), #overlap & leukemia
                            str(len(temp.loc[pd.notna(temp['TOPMed_CHIPVar'])])), #overlap & topmed
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')])), #overlap & somatic literature
                            #####
                            str(len(temp.query('Binomial_Prob < 0.01'))), #overlap & binom
                            str(len(temp.query('Binomial_Prob < 0.01 & LeukemiaGene == 1'))), #overlap & leukemia & binom
                            str(len(temp.loc[pd.isna(temp['TOPMed_CHIPVar'])].query('Binomial_Prob < 0.01'))), #overlap & topmed & binom
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Binomial_Prob < 0.01'))), #overlap & somatic literature & binomial
                            #####
                            str(len(temp.query('VAF > 0.1'))), #overlap & VAF
                            str(len(temp.query('VAF > 0.1 & LeukemiaGene == 1'))), #overlap & leukemia & VAF
                            str(len(temp.loc[pd.isna(temp['TOPMed_CHIPVar'])].query('VAF > 0.1'))), #overlap & topmed & VAF
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('VAF > 0.1'))), #overlap & somatic literature & VAF
                            #####
                            str(len(temp.query('Binomial_Prob < 0.01 & VAF > 0.1'))), #overlap & binom
                            str(len(temp.query('Binomial_Prob < 0.01 & VAF > 0.1 & LeukemiaGene == 1'))), #overlap & leukemia & binom & VAF
                            str(len(temp.loc[pd.isna(temp['TOPMed_CHIPVar'])].query('Binomial_Prob < 0.01 & VAF > 0.1'))), #overlap & topmed & binom & VAF
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Binomial_Prob < 0.01 & VAF > 0.1'))) #overlap & somatic literature & binomial & VAF
                           ]))
        fh.write('\n')
