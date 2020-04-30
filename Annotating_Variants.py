#----------------
#Maryam Zekavat
#Annotating CHUD variants in Hail (.py environment): 
conda activate hail
cd /Users/mzekavat/opt/anaconda3/envs/hail
hailctl dataproc start mz01 --master-machine-type n1-highmem-16 --worker-machine-type n1-highmem-16 --worker-boot-disk-size 200 --num-workers 10 --num-preemptible-workers 20 --master-boot-disk-size 100 --region us-east1 --zone us-east1-d --requester-pays-allow-all --vep GRCh38 --properties "spark:spark.driver.memory=90G,spark:spark.driver.maxResultSize=50G,spark:spark.kryoserializer.buffer.max=1G,spark:spark.task.maxFailures=20,spark:spark.driver.extraJavaOptions=-Xss4M,spark:spark.executor.extraJavaOptions=-Xss4M,spark:spark.speculation=true"
hailctl dataproc connect mz01 notebook --zone us-east1-d --region us-east1



import hail as hl
from pprint import pprint
from bokeh.io import output_notebook,show,save
from bokeh.layouts import gridplot
from bokeh.models import Span
import hail.expr.aggregators as agg
from bokeh.plotting import figure, output_file
import numpy as np
​
​
hl.init(default_reference='GRCh38',min_block_size=6)
​
​#Annotations: gsutil -m cp /medpop/esp2/mzekavat/CHIP/CHUD/data/variant_annot/somVariants.txt.bgz gs://maryam_lipids/UKBB_CHIP/somVariants.txt.bgz

kt = hl.import_table('gs://maryam_lipids/UKBB_CHIP/somVariants.txt.bgz', impute = True,min_partitions=2000,no_header = True) 
kt2 = kt.key_by(**hl.parse_variant(kt.f0)
kt2.describe()
kt2.write('gs://maryam_lipids/UKBB_CHIP/all_somatic_var_list.ht')

kt2=hl.read_table('gs://maryam_lipids/UKBB_CHIP/all_somatic_var_list.ht').repartition(1000)
kt2 = hl.vep(kt2, 'gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json')

consequence_in_severity_order = [
  "transcript_ablation"
, "splice_acceptor_variant"
, "splice_donor_variant"
, "stop_gained"
, "frameshift_variant"
, "stop_lost"
, "start_lost"
, "transcript_amplification"
, "inframe_insertion"
, "inframe_deletion"
, "missense_variant"
, "protein_altering_variant"
, "splice_region_variant"
, "incomplete_terminal_codon_variant"
, "stop_retained_variant"
, "synonymous_variant"
, "coding_sequence_variant"
, "mature_miRNA_variant"
, "5_prime_UTR_variant"
, "3_prime_UTR_variant"
, "non_coding_transcript_exon_variant"
, "intron_variant"
, "NMD_transcript_variant"
, "non_coding_transcript_variant"
, "upstream_gene_variant"
, "downstream_gene_variant"
, "TFBS_ablation"
, "TFBS_amplification"
, "TF_binding_site_variant"
, "regulatory_region_ablation"
, "regulatory_region_amplification"
, "feature_elongation"
, "regulatory_region_variant"
, "feature_truncation"
, "intergenic_variant"
]



mt4 = kt2.annotate(transcript_canonicals =kt2.vep.transcript_consequences.filter(lambda tc: tc.canonical == 1))

mt4 = mt4.annotate(
   all_transcript_terms = hl.set(hl.flatten(mt4.transcript_canonicals.map(lambda x: x.consequence_terms)))
)

mt4 = mt4.annotate(
    Consequence = hl.coalesce(  # coalesce means take first non-missing
        hl.literal(consequence_in_severity_order).filter(lambda cnsq:
            mt4.all_transcript_terms.contains(cnsq)
        ).head(),
        mt4.vep.most_severe_consequence
    )
)

mt4 = mt4.annotate(
	Gene = hl.if_else(mt4.transcript_canonicals.any(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)), \
				mt4.transcript_canonicals.find(lambda tc: (tc.canonical == 1) & (hl.set(tc.consequence_terms).contains(mt4.Consequence))).gene_symbol, \
			 mt4.vep.transcript_consequences.find(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)).gene_symbol), \

	hgvsp = hl.if_else(mt4.transcript_canonicals.any(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)), \
				mt4.transcript_canonicals.find(lambda tc: (tc.canonical == 1) & (hl.set(tc.consequence_terms).contains(mt4.Consequence))).hgvsp, \
			 mt4.vep.transcript_consequences.find(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)).hgvsp),

	hgvsc = hl.if_else(mt4.transcript_canonicals.any(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)), \
				mt4.transcript_canonicals.find(lambda tc: (tc.canonical == 1) & (hl.set(tc.consequence_terms).contains(mt4.Consequence))).hgvsc, \
			 mt4.vep.transcript_consequences.find(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)).hgvsc),
	
	SIFT = hl.if_else(mt4.transcript_canonicals.any(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)), \
				mt4.transcript_canonicals.find(lambda tc: (tc.canonical == 1) & (hl.set(tc.consequence_terms).contains(mt4.Consequence))).sift_prediction, \
			 mt4.vep.transcript_consequences.find(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)).sift_prediction),
	
	PolyPhen = hl.if_else(mt4.transcript_canonicals.any(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)), \
				mt4.transcript_canonicals.find(lambda tc: (tc.canonical == 1) & (hl.set(tc.consequence_terms).contains(mt4.Consequence))).polyphen_prediction, \
			 mt4.vep.transcript_consequences.find(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)).polyphen_prediction),

	LOF_LOFTEE = hl.if_else(mt4.transcript_canonicals.any(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)), \
                mt4.transcript_canonicals.find(lambda tc: (tc.canonical == 1) & (hl.set(tc.consequence_terms).contains(mt4.Consequence))).lof, \
             mt4.vep.transcript_consequences.find(lambda tc: hl.set(tc.consequence_terms).contains(mt4.Consequence)).lof),

	clin_sig = mt4.vep.colocated_variants.clin_sig,

	phenotype_or_disease = mt4.vep.colocated_variants.phenotype_or_disease,

	pubmed = mt4.vep.colocated_variants.pubmed,

	SOMATIC = mt4.vep.colocated_variants.somatic

	)

db = hl.experimental.DB()
mt4 = db.annotate_rows_db(mt4, "dbNSFP_variants")
#print(mt2.describe())
mt4 = mt4.annotate(metasvm = mt4.dbNSFP_variants.MetaSVM_pred)
mt4 = mt4.drop('dbNSFP_variants')

Leukemia_Genes = hl.import_table('gs://maryam_lipids/UKBB_CHIP/Leukemia_Genes.txt',no_header=False)
Leukemia_Genes = Leukemia_Genes.annotate(LeukemiaGene = 1).key_by('Gene_name')
mt4 = mt4.annotate(**Leukemia_Genes.index(mt4.Gene))

TOPMed_CHIPVar = hl.import_table('gs://maryam_lipids/UKBB_CHIP/topmed_CHIP_var.txt',no_header=False, impute=True)
TOPMed_CHIPVar = TOPMed_CHIPVar.annotate(CHROMv2 = ('chr'+TOPMed_CHIPVar.CHROM))
TOPMed_CHIPVar = TOPMed_CHIPVar.annotate(locus = hl.locus(TOPMed_CHIPVar.CHROMv2, TOPMed_CHIPVar.POS, reference_genome='GRCh38'), alleles = hl.array([TOPMed_CHIPVar.REF, TOPMed_CHIPVar.ALT]))
TOPMed_CHIPVar = TOPMed_CHIPVar.key_by(TOPMed_CHIPVar.locus, TOPMed_CHIPVar.alleles)

mt5 = mt4.annotate_rows(TOPMed_CHIPVar=TOPMed_CHIPVar[mt4.key])

#### Gnomad genomes annotations
gnomad_genomes = hl.read_table('gs://gnomad-public/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht') #.select('freq','popmax')
adj_nfe_index = gnomad_genomes.freq_index_dict['adj_nfe'].collect()[0]
adj_fin_index = gnomad_genomes.freq_index_dict['adj_fin'].collect()[0]
adj_oth_index = gnomad_genomes.freq_index_dict['adj_oth'].collect()[0]
adj_asj_index = gnomad_genomes.freq_index_dict['adj_asj'].collect()[0]
adj_afr_index = gnomad_genomes.freq_index_dict['adj_afr'].collect()[0]
adj_eas_index = gnomad_genomes.freq_index_dict['adj_eas'].collect()[0]
adj_sas_index = gnomad_genomes.freq_index_dict['adj_sas'].collect()[0]
adj_amr_index = gnomad_genomes.freq_index_dict['adj_amr'].collect()[0]
raw_index = gnomad_genomes.freq_index_dict['raw'].collect()[0]
adj_index = gnomad_genomes.freq_index_dict['adj'].collect()[0]
joined_genomes = gnomad_genomes[mt5.key]

mt6 = mt5.annotate(nfe_AF_g=joined_genomes.freq[adj_nfe_index].AF,fin_AF_g=joined_genomes.freq[adj_fin_index].AF,oth_AF_g=joined_genomes.freq[adj_oth_index].AF,asj_AF_g=joined_genomes.freq[adj_asj_index].AF,afr_AF_g=joined_genomes.freq[adj_afr_index].AF,eas_AF_g=joined_genomes.freq[adj_eas_index].AF,sas_AF_g=joined_genomes.freq[adj_sas_index].AF,amr_AF_g=joined_genomes.freq[adj_amr_index].AF,raw_AF_g=joined_genomes.freq[raw_index].AF,adj_AF_g=joined_genomes.freq[adj_index].AF)
#------------------------------------------------------------------------------------------------------------------------------------------------------------ #

####gnomad exomes annotations
gnomad_exomes = hl.read_table('gs://gnomad-public/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht') #.select('freq','popmax')
gnomad_nfe_index = gnomad_exomes.freq_index_dict['gnomad_nfe'].collect()[0]
gnomad_fin_index = gnomad_exomes.freq_index_dict['gnomad_fin'].collect()[0]
gnomad_oth_index = gnomad_exomes.freq_index_dict['gnomad_oth'].collect()[0] 
gnomad_asj_index = gnomad_exomes.freq_index_dict['gnomad_asj'].collect()[0]
gnomad_afr_index = gnomad_exomes.freq_index_dict['gnomad_afr'].collect()[0]
gnomad_eas_index = gnomad_exomes.freq_index_dict['gnomad_eas'].collect()[0]
gnomad_sas_index = gnomad_exomes.freq_index_dict['gnomad_sas'].collect()[0]
gnomad_amr_index = gnomad_exomes.freq_index_dict['gnomad_amr'].collect()[0]
raw_index = gnomad_exomes.freq_index_dict['gnomad_raw'].collect()[0]
joined_exomes = gnomad_exomes[mt6.key]

mt7 = mt6.annotate(nfe_AF_e=joined_exomes.freq[gnomad_nfe_index].AF,fin_AF_e=joined_exomes.freq[gnomad_fin_index].AF,oth_AF_e=joined_exomes.freq[gnomad_oth_index].AF,asj_AF_e=joined_exomes.freq[gnomad_asj_index].AF,afr_AF_e=joined_exomes.freq[gnomad_afr_index].AF,eas_AF_e=joined_exomes.freq[gnomad_eas_index].AF,sas_AF_e=joined_exomes.freq[gnomad_sas_index].AF,amr_AF_e=joined_exomes.freq[gnomad_amr_index].AF,raw_AF_e=joined_exomes.freq[raw_index].AF)

LOF_Conseq = ["frameshift_variant", "transcript_ablation" , "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "start_lost"]
mt8 = mt7.filter((((((hl.literal(LOF_Conseq).contains(mt7.Consequence)) | ((mt7.metasvm).contains('D')) | (mt7.LOF_LOFTEE.contains("HC") ))) & \
                                                ((mt7.raw_AF_e <= 0.01) | (~(hl.is_defined(mt7.raw_AF_e)))) & \
                                                ((mt7.afr_AF_e < 0.01)  | (~(hl.is_defined(mt7.afr_AF_e)))) & \
                                                ((mt7.amr_AF_e < 0.01) | (~(hl.is_defined(mt7.amr_AF_e)))) & \
                                                ((mt7.asj_AF_e < 0.01) | (~(hl.is_defined(mt7.asj_AF_e)))) & \
                                                ((mt7.eas_AF_e < 0.01) | (~(hl.is_defined(mt7.eas_AF_e)))) & \
                                                ((mt7.fin_AF_e < 0.01) | (~(hl.is_defined(mt7.fin_AF_e)))) & \
                                                ((mt7.nfe_AF_e < 0.01) | (~(hl.is_defined(mt7.nfe_AF_e)))) & \
                                                ((mt7.sas_AF_e < 0.01) | (~(hl.is_defined(mt7.sas_AF_e)))) & \
                                                ((mt7.oth_AF_e < 0.01) | (~(hl.is_defined(mt7.oth_AF_e)))) & \
                                                ((mt7.raw_AF_g <= 0.01) | (~(hl.is_defined(mt7.raw_AF_g)))) & \
                                                ((mt7.afr_AF_g < 0.01)  | (~(hl.is_defined(mt7.afr_AF_g)))) & \
                                                ((mt7.amr_AF_g < 0.01) | (~(hl.is_defined(mt7.amr_AF_g)))) & \
                                                ((mt7.asj_AF_g < 0.01) | (~(hl.is_defined(mt7.asj_AF_g)))) & \
                                                ((mt7.eas_AF_g < 0.01) | (~(hl.is_defined(mt7.eas_AF_g )))) & \
                                                ((mt7.fin_AF_g < 0.01) | (~(hl.is_defined(mt7.fin_AF_g)))) & \
                                                ((mt7.nfe_AF_g < 0.01) | (~(hl.is_defined(mt7.nfe_AF_g)))) & \
                                                ((mt7.oth_AF_g < 0.01) | (~(hl.is_defined(mt7.oth_AF_g))))) | (hl.is_defined(mt7.TOPMed_CHIPVar.Gene))), keep = True) 
mt8= mt8.annotate(Exon = hl.if_else(mt8.transcript_canonicals.any(lambda tc: hl.set(tc.consequence_terms).contains(mt8.Consequence)), \
				mt8.transcript_canonicals.find(lambda tc: (tc.canonical == 1) & (hl.set(tc.consequence_terms).contains(mt8.Consequence))).exon, \
			 mt8.vep.transcript_consequences.find(lambda tc: hl.set(tc.consequence_terms).contains(mt8.Consequence)).exon))

mt8 = mt8.drop('vep', 'transcript_canonicals')
mt8.write('gs://maryam_lipids/UKBB_CHIP/somVariants.filtered_rare_disruptive_LOF.annotated.ht') #1409376 variants among this list

#now count the # variants that are among the Leukemia Genes
mt8=hl.read_table('gs://maryam_lipids/UKBB_CHIP/somVariants.filtered_rare_disruptive_LOF.annotated.ht')

print(mt8.filter((mt8.LeukemiaGene == 1)).count()) #12348 variants
#now count the # variants that were seen in TOPMed 
print(mt8.filter((hl.is_defined(mt8.TOPMed_CHIPVar.Gene))).count()) #613 variants

​
