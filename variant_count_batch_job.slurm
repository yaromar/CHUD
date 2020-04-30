#!/bin/bash
#SBATCH --partition=general
#SBATCH --mem-per-cpu=10G
#SBATCH --time=6:00:00

module load miniconda/4.8.2
source activate py37_dev


BATCHNUM=1
BATCHTOT=10

python variantCountBatch.py /ysm-gpfs/project/smz25/CHUD/UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz /ysm-gpfs/project/smz25/CHUD/varcounts$BATCHNUM.txt /ysm-gpfs/project/smz25/CHUD/QC_v1/Filtered_SomaticCalls_v1/ $BATCHNUM/$BATCHTOT

# python variantCountBatch.py /ysm-gpfs/project/smz25/CHUD/UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz /ysm-gpfs/project/smz25/CHUD/varcounts1.txt /ysm-gpfs/project/smz25/CHUD/QC_v1/Filtered_SomaticCalls_v1/ 1/10

