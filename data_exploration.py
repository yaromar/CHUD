import pandas as pd

testfile = "UKBB_CHIP-somVariants.filtered_rare_disruptive_LOF.annotated.bgz"

test_data = pd.read_csv(testfile, delimiter='\t', compression='gzip')


