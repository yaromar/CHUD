
# Success. You can now start the database server using:

#     postgres -D chuddb
# or
#     pg_ctl -D chuddb -l logfile start



# psql -h localhost -p 5432 chuddb

# COPY Variant FROM '/ysm-gpfs/project/smz25/CHUD/UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.txt' DELIMITER E'\t';

import pandas as pd

# var_annotations = "UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz"
# annotations = pd.read_csv(var_annotations, delimiter='\t', compression='gzip')
# print(annotations.columns)
# print(annotations.shape)
# print(annotations.dtypes)
# print(annotations.head())

# print("\n\n")


# variant_file = "2657375_23173_0_0-filtered.vcf.gz.1343.GTs.bgz"
# variants = pd.read_csv(variant_file, delimiter='\t', compression='gzip')
# print(variants.columns)
# print(variants.shape)
# print(variants.dtypes)
# print(variants.head())

# print("\n\n")

pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 200)

phenotype_file = "ukbb_PhenoFile.ForClassProj.QCed.txt.gz"
phenotypes = pd.read_csv(phenotype_file, delimiter='\t', compression='gzip')
print(phenotypes.columns)
print(phenotypes.shape)
print(phenotypes.dtypes)
print(phenotypes.head())

