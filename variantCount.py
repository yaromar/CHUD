import pandas as pd
from os import listdir
import numpy as np

#the annotation file should be in the same directory as the script
annotations = pd.read_csv('projects_mzekavat_CHIP_somVariants.filtered_rare_disruptive_LOF.annotated.plusLCR.plusAnnovar.bgz', delimiter='\t', compression='gzip')[['f0', 'SOMATIC', 'TOPMed_CHIP_Var', 'UKBB_CHIP_wl_Var', 'UKBB_DNMT3A_otherMis_Var', 'Known_400Leukemia_Gene', 'Known_74CHIP_Gene', 'Hematopoietic_COSMIC', 'LCR']].set_index('f0')

annotationsLCR = annotations.fillna(0).query('LCR == 0')
annotationsLCR.at['chr2:197401887:C:T', 'UKBB_CHIP_wl_Var'] = 1
annotations = annotations.fillna({'SOMATIC': '0'})
annotationsLCR = annotations.fillna(0)

#the FOLDER with the somatic calls files should be in the current directory
#the files should be compressed and have .bgz extension
#there should be no other files
with open('variantCount_output.txt', 'w') as fh:
    fh.write('\t'.join(['sample_id',
                        'overlap & UKBB CHIP',
                        'overlap & (topmed OR UKBB CHIP)',
                        'overlap & (topmed OR ukbb chip OR ukbb dnmt3a)',

                        'mutectVar',
                        'overlap',
                        'overlap & leukemia gene',
                        'overlap & 74chip gene',

                        'mutectVar & binom',
                        'overlap & binom',
                        'overlap & leukemia gene & binom',
                        'overlap & 74chip gene & binom',

                        'mutectVar & vaf',
                        'overlap & vaf',
                        'overlap & leukemia gene & vaf',
                        'overlap & 74chip gene & vaf',

                        'mutectVar & vaf & binom',
                        'overlap & vaf & binom',
                        'overlap & leukemia gene & vaf & binom',
                        'overlap & 74chip gene & vaf & binom',



                        'overlap & somatic literature',
                        'overlap & leukemia gene & somatic literature',
                        'overlap & 74chip gene & somatic literature',

                        'overlap & somatic literature & binom',
                        'overlap & leukemia gene & somatic literature & binom',
                        'overlap & 74chip gene & somatic literature & binom',

                        'overlap & somatic literature & vaf',
                        'overlap & leukemia gene & somatic literature & vaf',
                        'overlap & 74chip gene & somatic literature & vaf',

                        'overlap & somatic literature & vaf & binom',
                        'overlap & leukemia gene & somatic literature & vaf & binom',
                        'overlap & 74chip gene & somatic literature & vaf & binom',



                        'overlap & Hematopoietic_COSMIC',
                        'overlap & leukemia gene & Hematopoietic_COSMIC',
                        'overlap & 74chip gene & Hematopoietic_COSMIC',

                        'overlap & Hematopoietic_COSMIC & binom',
                        'overlap & leukemia gene & Hematopoietic_COSMIC & binom',
                        'overlap & 74chip gene & Hematopoietic_COSMIC & binom',

                        'overlap & Hematopoietic_COSMIC & vaf',
                        'overlap & leukemia gene & Hematopoietic_COSMIC & vaf',
                        'overlap & 74chip gene & Hematopoietic_COSMIC & vaf',

                        'overlap & Hematopoietic_COSMIC & vaf & binom',
                        'overlap & leukemia gene & Hematopoietic_COSMICe & vaf & binom',
                        'overlap & 74chip gene & Hematopoietic_COSMIC & vaf & binom',



                        'yaro'
						]))
    fh.write('\n')

    for file in listdir('Sparse_QCplusFilterMutectPass'):
        sample = pd.read_csv('Sparse_QCplusFilterMutectPass/'+file, delimiter='\t', compression='gzip', usecols=["v", "Binomial_Prob", "VAF"]).set_index('v')
        temp = pd.concat([annotationsLCR, sample], axis=1, join='inner')

        sampleID, mutectVar = [file.split('.', 4)[i] for i in (0, 3)]
        sampleID = sampleID.split('.', 4)[0].split('-')[0]

        fh.write('\t'.join([sampleID,
                            str(len(temp.query('UKBB_CHIP_wl_Var == 1'))), #overlap & UKBB CHIP
                            str(len(temp.query('UKBB_CHIP_wl_Var == 1 | TOPMed_CHIP_Var == 1'))), #overlap & (topmed OR UKBB CHIP)
                            str(len(temp.query('UKBB_CHIP_wl_Var == 1 | TOPMed_CHIP_Var == 1 | UKBB_DNMT3A_otherMis_Var == 1'))), #overlap & (topmed OR ukbb chip OR ukbb dnmt3a)

                            str(len(sample)), #mutectVar
                            str(len(temp)), #overlap
                            str(len(temp.query('Known_400Leukemia_Gene == 1'))), #overlap & leukemia gene
                            str(len(temp.query('Known_74CHIP_Gene == 1'))), #overlap & 74chip gene

                            str(len(sample.query('Binomial_Prob < 0.01'))), #mutectVar & binom
                            str(len(temp.query('Binomial_Prob < 0.01'))), #overlap & binom
                            str(len(temp.query('Known_400Leukemia_Gene == 1 & Binomial_Prob < 0.01'))), #overlap & leukemia gene & binom
                            str(len(temp.query('Known_74CHIP_Gene == 1 & Binomial_Prob < 0.01'))), #overlap & 74chip gene & binom

                            str(len(sample.query('VAF > 0.1'))), #mutectVar & vaf
                            str(len(temp.query('VAF > 0.1'))), #overlap & vaf
                            str(len(temp.query('Known_400Leukemia_Gene == 1 & VAF > 0.1'))), #overlap & leukemia gene & vaf
                            str(len(temp.query('Known_74CHIP_Gene == 1 & VAF > 0.1'))), #overlap & 74chip gene & vaf

                            str(len(sample.query('VAF > 0.1 & Binomial_Prob < 0.01'))), #mutectVar & vaf & binom
                            str(len(temp.query('VAF > 0.1 & Binomial_Prob < 0.01'))), #overlap & vaf & binom
                            str(len(temp.query('Known_400Leukemia_Gene == 1 & VAF > 0.1 & Binomial_Prob < 0.01'))), #overlap & leukemia gene & vaf & binom
                            str(len(temp.query('Known_74CHIP_Gene == 1 & VAF > 0.1 & Binomial_Prob < 0.01'))), #overlap & 74chip gene & vaf & binom



                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')])), #overlap & somatic literature
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_400Leukemia_Gene == 1'))), #overlap & leukemia gene & somatic literature
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_74CHIP_Gene == 1'))), #overlap & 74chip gene & somatic literature

                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Binomial_Prob < 0.01'))), #overlap & somatic literature & binom
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_400Leukemia_Gene == 1 & Binomial_Prob < 0.01'))), #overlap & leukemia gene & somatic literature & binom
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_74CHIP_Gene == 1 & Binomial_Prob < 0.01'))), #overlap & 74chip gene & somatic literature & binom

                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('VAF > 0.1'))), #overlap & somatic literature & VAF
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_400Leukemia_Gene == 1 & VAF > 0.1'))), #overlap & leukemia gene & somatic literature & VAF
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_74CHIP_Gene == 1 & VAF > 0.1'))), #overlap & 74chip gene & somatic literature & VAF

                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('VAF > 0.1 & Binomial_Prob < 0.01'))), #overlap & somatic literature & VAF & binom
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_400Leukemia_Gene == 1 & VAF > 0.1 & Binomial_Prob < 0.01'))), #overlap & leukemia gene & somatic literature & VAF & binom
                            str(len(temp.loc[temp['SOMATIC'].str.contains('1')].query('Known_74CHIP_Gene == 1 & VAF > 0.1 & Binomial_Prob < 0.01'))), #overlap & 74chip gene & somatic literature & VAF & binom



                            str(len(temp)), #overlap & Hematopoietic_COSMIC & Hematopoietic_COSMIC
                            str(len(temp.query('Known_400Leukemia_Gene == 1 & Hematopoietic_COSMIC == 1'))), #overlap & leukemia gene & Hematopoietic_COSMIC
                            str(len(temp.query('Known_74CHIP_Gene == 1 & Hematopoietic_COSMIC == 1'))), #overlap & 74chip gene & Hematopoietic_COSMIC

                            str(len(temp.query('Binomial_Prob < 0.01 & Hematopoietic_COSMIC == 1'))), #overlap & binom & Hematopoietic_COSMIC
                            str(len(temp.query('Known_400Leukemia_Gene == 1 & Binomial_Prob < 0.01 & Hematopoietic_COSMIC == 1'))), #overlap & leukemia gene & binom & Hematopoietic_COSMIC
                            str(len(temp.query('Known_74CHIP_Gene == 1 & Binomial_Prob < 0.01 & Hematopoietic_COSMIC == 1'))), #overlap & 74chip gene & binom & Hematopoietic_COSMIC

                            str(len(temp.query('VAF > 0.1 & Hematopoietic_COSMIC == 1'))), #overlap & vaf & Hematopoietic_COSMIC
                            str(len(temp.query('Known_400Leukemia_Gene == 1 & VAF > 0.1 & Hematopoietic_COSMIC == 1'))), #overlap & leukemia gene & vaf & Hematopoietic_COSMIC
                            str(len(temp.query('Known_74CHIP_Gene == 1 & VAF > 0.1 & Hematopoietic_COSMIC == 1'))), #overlap & 74chip gene & vaf & Hematopoietic_COSMIC

                            str(len(temp.query('VAF > 0.1 & Binomial_Prob < 0.01 & Hematopoietic_COSMIC == 1'))), #overlap & vaf & binom & Hematopoietic_COSMIC
                            str(len(temp.query('Known_400Leukemia_Gene == 1 & VAF > 0.1 & Binomial_Prob < 0.01 & Hematopoietic_COSMIC == 1'))), #overlap & leukemia gene & vaf & binom & Hematopoietic_COSMIC
                            str(len(temp.query('Known_74CHIP_Gene == 1 & VAF > 0.1 & Binomial_Prob < 0.01 & Hematopoietic_COSMIC == 1'))), #overlap & 74chip gene & vaf & binom & Hematopoietic_COSMIC



                            ####YARO:
                            str(len(sample.query('(VAF > 0.02 & VAF < 0.46) | (VAF > 0.54 & VAF < 0.98)'))), #potential clonal expansion




                           ]))
        fh.write('\n')
