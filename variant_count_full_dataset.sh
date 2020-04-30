#!/bin/bash
# This is a comment!


BATCHTOT=10

for BATCHNUM in $(seq 1 $BATCHTOT)
# for BATCHNUM in 1 2 3 4 5 6 7 8 9 10
do
   sbatch --job-name=soc_varcount$BATCHNUM.$BATCHTOT.run --output=soc_varcount$BATCHNUM.$BATCHTOT.out --export=ALL,BATCHNUM=$BATCHNUM,BATCHTOT=$BATCHTOT variant_count_batch_job.slurm
   sleep 2
done

# dsq --job-file variantcountbatchjob.txt --job-name dsq-socvarcount --mem-per-cpu 10g -t 8:00:00 --mail-type ALL 

# BATCHTOT=10
# for BATCHNUM in $(seq 1 $BATCHTOT)
# do
# 	echo "module load miniconda; source activate py37_dev; python variantCountBatch.py /ysm-gpfs/project/smz25/CHUD/UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz /ysm-gpfs/project/smz25/CHUD/varcounts$BATCHNUM.txt /ysm-gpfs/project/smz25/CHUD/QC_v1/Filtered_SomaticCalls_v1/ $BATCHNUM/$BATCHTOT"
# done


