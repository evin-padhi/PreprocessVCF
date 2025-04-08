version 1.0

workflow PreprocessVCF {
    input {
        String vds_path
        File ancestry_labels
    }

    call SplitVDS {
        input:
            vds_path = vds_path
    }    

}

task SplitVDS {
    input {
        String vds_path

        String docker = "hailgenetics/hail:0.2.134-py3.11"
        Int cpu = 4
        Int memory_gb = 15
        Int disk_size = 150
    }

    command <<<
set -euo pipefail

# Write Python script to a file
cat << 'EOF' > split_vds.py
import hail as hl
import sys
import re

def make_dense_mt(vds: hl.vds.VariantDataset, is_filtering_FT: bool, is_keep_as_vqsr_fields: bool,
                  max_alt_alleles: Optional[int], fields_to_drop: str) -> hl.MatrixTable:
    """
    Convert a Variant Dataset (VDS) to a densified MatrixTable (MT) with additional annotations, and drop some
    fields if fields exist in the VDS.

    Parameters:
    vds (hl.vds.VariantDataset): The input Variant Dataset.
    is_filtering_FT (bool): Flag to filter GTs based on the FT field. If "FT" is "Fail", then set
                        corresponding "GT" to missing.
    is_keep_as_vqsr_fields (bool): Flag to keep VQSR fields.
    max_alt_alleles (Optional[int]): Maximum number of alternate alleles for filtering.
    fields_to_drop (str): Specific fields to drop from the VDS while converting to a dense MT, e.g., "as_vqsr, LAD, tranche_data,
                          truth_sensitivity_snp_threshold, truth_sensitivity_indel_threshold,
                          snp_vqslod_threshold, indel_vqslod_threshold", delimiter can be , _ or space.

    Returns:
    hl.MatrixTable: The resulting densified and annotated MatrixTable.
    """

    #vd_gt = to_dense_mt(vds)
    vd_gt = vds.variant_data

    if max_alt_alleles is not None:
        vd_gt = vd_gt.filter_rows(hl.len(vd_gt.alleles) < max_alt_alleles)

    vd_gt = vd_gt.annotate_entries(AD=hl.vds.local_to_global(vd_gt.LAD, vd_gt.LA, n_alleles=hl.len(vd_gt.alleles), fill_value=0, number='R'))
    vd_gt = vd_gt.transmute_entries(GT=hl.vds.lgt_to_gt(vd_gt.LGT, vd_gt.LA))

    if 'FT' in vd_gt.entry:
        vd_gt = vd_gt.transmute_entries(FT=hl.if_else(vd_gt.FT, "PASS", "FAIL"))

    if 'gvcf_info' in vd_gt.entry:
        vd_gt = vd_gt.drop('gvcf_info')

    d_callset = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, vd_gt))
    #d_callset = vd_gt
    if is_filtering_FT and "FT" in d_callset.entry:
        d_callset = d_callset.annotate_entries(GT=hl.or_missing((~hl.is_defined(d_callset.FT)) | (d_callset.FT.contains("PASS")), d_callset.GT))

    d_callset = hl.variant_qc(d_callset)
    d_callset = d_callset.annotate_rows(info=hl.struct(AC=d_callset.variant_qc.AC[1:], AF=d_callset.variant_qc.AF[1:], AN=d_callset.variant_qc.AN, homozygote_count=d_callset.variant_qc.homozygote_count))

    if is_keep_as_vqsr_fields and ('as_vqsr' in d_callset.row):
        d_callset = d_callset.annotate_rows(info=d_callset.info.annotate(AS_VQSLOD=d_callset.as_vqsr.values().vqslod, AS_YNG=d_callset.as_vqsr.values().yng_status))

    d_callset = d_callset.annotate_rows(variant_qc=d_callset.variant_qc.drop("AC", "AF", "AN", "homozygote_count"))
    fields_to_drop_list = re.split(r'[,\s]+', fields_to_drop)
    d_callset = d_callset.drop(*(f for f in fields_to_drop_list if f in d_callset.entry or f in d_callset.row or f in d_callset.col or f in d_callset.globals))

    return d_callset

def calculate_info(mt):
    mt = mt.drop("variant_qc", "info")
    mt = hl.variant_qc(mt)

    mt = mt.annotate_rows(info = hl.struct(AC=mt.variant_qc.AC[1:],
                                       AF=mt.variant_qc.AF[1:],
                                       AN=mt.variant_qc.AN,
                                       homozygote_count=mt.variant_qc.homozygote_count))
    
    # Drop duplicate nested fields that are already in INFO field rendered by call_stats()
   # mt = mt.annotate_rows(variant_qc = mt.variant_qc.drop("AC", "AF", "AN", "homozygote_count"))
    return mt

def split_mt(mt):
    # Filter rows that have a non-missing value.  
    #  Note that this will have an issue w/ "PASS"
    mt = mt.filter_rows(hl.is_missing(mt.filters) | (hl.len(mt['filters']) == 0), keep=True)
    split = hl.split_multi_hts(mt)
    split = calculate_info(split)
    return split

hl.init()

vds_path = sys.argv[1]
print(f"VDS path provided: {vds_path}")

if not vds_path:
    raise ValueError("VDS path argument is empty!")

vds = hl.vds.read_vds(vds_path)

chromosomes = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
vds_chromosomes = {chr: hl.vds.filter_chromosomes(vds, keep=chr) for chr in chromosomes}
mt_chromosomes = {chr: make_dense_mt(vds_chromosomes[chr]) for chr in chromosomes}

for chr, mt in mt_chromosomes.items():
    hl.export_vcf(mt, f'{chr}.vcf.bgz', tabix=True)

hl.stop()
EOF

# Run the Python script with the argument
python3 split_vds.py "~{vds_path}"
>>>

    runtime {
        docker: docker
        memory: memory_gb + ' GB'
        cpu: cpu
        disk: "local-disk ~{disk_size} HDD"
        }

    output {
        Array[File] vcfs = glob("*.vcf.bgz")
    }
}


# Task to process each chromosome
task ProcessChromosome {
    input {
        File vcf_file
        File vcf_file_index
        File ancestry_labels
        String min_allele_count
        String chromosome

        String docker = "quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
        Int cpu = 1
        Int memory_gb = ceil(size(vcf_file, "GB") + 8)
        Int disk_size = ceil(size(vcf_file, "GB") + 20)
    }

    command {
        # Use bcftools to filter the VCF for the specific chromosome
        bcftools view -r ${chromosome} ${vcf_file} -Ou \
        | bcftools norm -m + - \
        | bcftools view -m2 -M2 -c${min_allele_count} - \
        | bcftools +fill-tags -  -- -t AN,AC,AF,MAF,AC_hom,AC_het,HWE,F_MISSING -S ${ancestry_labels} \
        | bcftools view -e "F_MISSING > 0.05" -W -Oz -o ${chromosome}.filtered.vcf.gz
    }
    
    runtime {
        docker: docker
        memory: memory_gb + ' GB'
        cpu: cpu
        disks: disk_size + ' GB'
    }
    
    output {
        File vcf_output = "${chromosome}.filtered.vcf.gz"
        File vcf_output_index = "${chromosome}.filtered.vcf.gz.csi"
    }
}

# Task to gather results (optional, depending on your needs)
task GatherResults {
    input {
        Array[File] chromosome_results
        Array[File] indices

        String docker = "quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
        Int cpu = 4
        Int memory_gb = 64
        Int disk_size = ceil(size(chromosome_results, "GB") + 20)
    }

    command {
        # Combine the filtered VCFs into one file (optional)
        bcftools concat --threads ~{cpu - 1} -Oz -o combined.vcf.gz ~{sep=' ' chromosome_results}
    }
    
    runtime {
        docker: docker
        memory: memory_gb + ' GB'
        cpu: cpu
        disk: disk_size + ' GB'
    }
    
    output {
        File combined_vcf = "combined.vcf.gz"
    }
}










