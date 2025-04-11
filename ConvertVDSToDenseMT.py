import hail as hl
import sys

def make_dense_mt(vds, filter_FT, keep_as_vqsr, max_alt, fields_drop):
    vd_gt = vds.variant_data
    if max_alt is not None:
        vd_gt = vd_gt.filter_rows(hl.len(vd_gt.alleles) < max_alt)
    vd_gt = vd_gt.annotate_entries(AD=hl.vds.local_to_global(vd_gt.LAD, vd_gt.LA, n_alleles=hl.len(vd_gt.alleles), fill_value=0, number='R'))
    vd_gt = vd_gt.transmute_entries(GT=hl.vds.lgt_to_gt(vd_gt.LGT, vd_gt.LA))
    if 'FT' in vd_gt.entry:
        vd_gt = vd_gt.transmute_entries(FT=hl.if_else(vd_gt.FT, "PASS", "FAIL"))
    if 'gvcf_info' in vd_gt.entry:
        vd_gt = vd_gt.drop('gvcf_info')
    d_callset = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, vd_gt))
    if filter_FT and "FT" in d_callset.entry:
        d_callset = d_callset.annotate_entries(GT=hl.or_missing((~hl.is_defined(d_callset.FT)) | (d_callset.FT.contains("PASS")), d_callset.GT))
    d_callset = hl.variant_qc(d_callset)
    d_callset = d_callset.annotate_rows(info=hl.struct(AC=d_callset.variant_qc.AC[1:], AF=d_callset.variant_qc.AF[1:], AN=d_callset.variant_qc.AN, homozygote_count=d_callset.variant_qc.homozygote_count))
    fields_drop_list = fields_drop.replace(",", " ").split()
    d_callset = d_callset.drop(*fields_drop_list)
    return d_callset


vds_url = sys.argv[1]
rnaseq_samples_tsv = sys.argv[2]
output_bucket = sys.argv[3]
identifier = sys.argv[4]

hl.init(default_reference='GRCh38')
vds = hl.vds.read_vds("${vds_url}")
fields_to_drop = "as_vqsr,LAD,LGT,LA,tranche_data,truth_sensitivity_snp_threshold,truth_sensitivity_indel_threshold,snp_vqslod_threshold,indel_vqslod_threshold"
samples = hl.import_table("${rnaseq_samples_tsv}", key="research_id")
vds_subset = hl.vds.filter_samples(vds, samples, keep=True, remove_dead_alleles=True)
mt = make_dense_mt(vds_subset, True, False, 100, fields_to_drop)
mt.write("${output_bucket}/${identifier}.mt", overwrite=True)
