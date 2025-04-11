import hail as hl
import sys

input_mt_path=sys.argv[1]
output_bucket=sys.argv[2]
identifier=sys.arv[3]

hl.init(default_reference='GRCh38')
mt = hl.read_matrix_table(input_mt_path)
hl.export_vcf(mt, f"{output_bucket}/{identifier}.vcf.bgz")
