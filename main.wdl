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

python3 <<CODE
import hail as hl

hl.init()

vds_path = os.getenv("VDS_PATH")
vds = hl.vds.read_vds(vds_path)

chromosomes = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
vds_chromosomes = {chr: hl.vds.filter_chromosomes(vds, keep=chr) for chr in chromosomes}
mt_chromosomes = {chr: hl.vds.to_dense_mt(vds_chromosomes[chr]) for chr in chromosomes}

for chr, mt in mt_chromosomes.items():
    hl.export_vcf(mt, f'{chr}.vcf.bgz', tabix=True)

hl.stop()
CODE
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










