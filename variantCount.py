import pandas as pd
from os import listdir
import numpy as np


annotations = pd.read_csv('UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz', delimiter='\t', compression='gzip')[['f0', 'SOMATIC', 'LeukemiaGene', 'TOPMed_CHIPVar']].set_index('f0')
annotations = annotations.fillna({'SOMATIC': '0'})

with open('variantCount_output.txt', 'w') as fh:
    fh.write('\t'.join(['sample_id', 'num_FilterMutect', 'annotated_overlap', 'SOMATIC', 'LeukemiaGene', 'TOPMed_CHIPVar', '\n']))
    
    for file in listdir('Filtered_SomaticCalls_v1'):
        sample = pd.read_csv('Filtered_SomaticCalls_v1/'+file, delimiter='\t', compression='gzip', usecols=[0]).set_index('v')
        temp = pd.concat([annotations, sample], axis=1, join='inner')
        sampleID, mutectVar = [file.split('.', 4)[i] for i in (0, 3)]
        sampleID = sampleID.split('.', 4)[0].split('-')[0]

        fh.write('\t'.join([sampleID, mutectVar, str(len(pd.concat([annotations, sample], axis=1, join='inner').index)), str(len(temp.loc[temp['LeukemiaGene'] == 1].index)), str(len(temp.loc[pd.notna(temp['TOPMed_CHIPVar'])].index)), str(len(temp.loc[temp['SOMATIC'].str.contains('1')])), '\n']))
