import hail as hl
import sys

input_mt_path=sys.argv[1]
output_bucket=sys.argv[2]
identifier=sys.argv[3]

hl.init(default_reference='GRCh38')
mt = hl.read_matrix_table(input_mt_path)
mt = mt.filter_rows(hl.is_missing(mt.filters) | (hl.len(mt.filters) == 0))
mt_split = hl.split_multi_hts(mt)
mt_split = hl.variant_qc(mt_split)
mt_split = mt_split.annotate_rows(info=hl.struct(AC=mt_split.variant_qc.AC[1:], AF=mt_split.variant_qc.AF[1:], AN=mt_split.variant_qc.AN, homozygote_count=mt_split.variant_qc.homozygote_count))
mt_split.write(f"{output_bucket}/{identifier}_split.mt", overwrite=True)
