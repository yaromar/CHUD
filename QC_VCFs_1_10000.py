
#conda activate hail
#cd /Users/mzekavat/opt/anaconda3/envs/hail
#hailctl dataproc start mz02 --master-machine-type n1-highmem-16 --worker-machine-type n1-highmem-16 --worker-boot-disk-size 200 --num-workers 3 --num-preemptible-workers 3 --master-boot-disk-size 100 --region us-east1 --zone us-east1-d --requester-pays-allow-all --properties "spark:spark.driver.memory=90G,spark:spark.driver.maxResultSize=50G,spark:spark.kryoserializer.buffer.max=1G,spark:spark.task.maxFailures=20,spark:spark.driver.extraJavaOptions=-Xss4M,spark:spark.executor.extraJavaOptions=-Xss4M,spark:spark.speculation=true"
#hailctl dataproc connect mz02 notebook --zone us-east1-d --region us-east1
#hailctl dataproc submit --zone us-east1-d --region us-east1 mz02 ~/Documents/Broad_2015_17/Python_Scripts_Hail/CHIP/Merge_SomaticVCFS_15000_30000.py
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference = "GRCh38")
import numpy as np
import pandas as pd
from collections import Counter
from math import log, isnan
from pprint import pprint
import time
from bokeh.io import show, output_notebook
from bokeh.layouts import gridplot
output_notebook()

recoding_dict = {f"{i + 1}":  f"chr{i + 1}" for i in range(22)}
recoding_dict['X'] = 'chrX'
recoding_dict['Y'] = 'chrY'

files = hl.import_table('gs://maryam_lipids/UKBB_CHIP/filenames.txt',impute=True,no_header=True)
files_list = [row['f0'] for row in files.select(files.f0).collect()]

for num in range(1,10000):
	print(num)
	filenamev2 = files_list[num].strip()
	mt=hl.import_vcf(filenamev2,filter='\d/\d/\d', skip_invalid_loci=True,force_bgz=True,reference_genome='GRCh38',contig_recoding=recoding_dict)
	### filter to pass variants and split_multi
	mt2 = mt.filter_rows(mt.filters.size() > 0,keep = False)
	mt2 = hl.split_multi_hts(mt2)

	#variant read counts of 3
	#at least one read in both forward and reverse orientations
	#remove  monomorphic variants
	mt3 = mt2.filter_entries(((mt2.AD[1]<2) | (mt2.F1R2[1]==0) | (mt2.F2R1[1]==0)) ,keep = False)
	mt3 = hl.variant_qc(mt3)
	mt3 = mt3.filter_rows( (mt3.variant_qc.AF[1] > 0) & (mt3.variant_qc.AF[1] < 1) ,keep = True)

	mt4 = mt3.annotate_rows(v = hl.variant_str(mt3.locus, mt3.alleles),\
						NumAltAlleles = hl.agg.max(mt3.GT.n_alt_alleles()), \
						VAF =hl.agg.explode(lambda x: hl.agg.mean(x), mt3.AF),\
						TLOD =mt3.info.TLOD[0], \
						GERMQ = mt3.info.GERMQ, \
						STR=mt3.info.STR,\
						AD_alt=hl.agg.mean(mt3.AD[1]),\
						AD_ref=hl.agg.mean(mt3.AD[0]))

	mt4 = mt4.annotate_entries(Binomial_Prob = hl.binom_test(mt4.AD[1],mt4.DP,0.5,'greater'))
	mt4 = mt4.key_rows_by("v")
	mt4 = mt4.drop('locus', 'alleles', 'qual', 'filters','variant_qc', 'GQ', 'PGT', 'PID', 'PL', 'PS', 'info', 'rsid', 'a_index', 'was_split')
	filt2=mt4.count_rows()
	mt4.entries().export(filenamev2+"."+str(filt2)+".GTs.bgz")
	del(mt)
	del(mt2)
	del(mt3)
	del(mt4)
