version 1.0

workflow WriteVCFWorkflow {
    input {
        String matrix_table
        String samples_table
        String ancestry_table
        String ancestry
        String chr
        Int MinimumAC_inclusive
        String output_path
    }

    call WriteVCFTask {
        input:
            matrix_table = matrix_table,
            samples_table = samples_table,
            ancestry_table = ancestry_table,
            ancestry = ancestry,
            chr = chr,
            MinimumAC_inclusive = MinimumAC_inclusive,
            output_path = output_path
    }

    call plink2 {
        input:
            vcf_file = WriteVCFTask.output_vcf,
            prefix = prefix,
            new_id_max_allele_len = new_id_max_allele_len
    }

    output {
        File output_vcf = WriteVCFTask.output_vcf
        Directory plink_outputs = plink2.plink_outputs
    }
}

task WriteVCFTask {
    input {
        String matrix_table
        String samples_table
        String ancestry_table
        String ancestry
        String chr
        Int MinimumAC_inclusive
        String output_path
    }

    command <<<
        set -e

        curl -O https://raw.githubusercontent.com/jonnguye/PreprocessVCF/NotebookToWDL/write_vcf.py

        python3 write_vcf.py \
            --matrix_table "~{matrix_table}" \
            --samples_table "~{samples_table}" \
            --ancestry_table "~{ancestry_table}" \
            --ancestry "~{ancestry}" \
            --chr "~{chr}" \
            --MinimumAC_inclusive "~{MinimumAC_inclusive}" \
            --output_path "~{output_path}"
    >>>

    runtime {
        docker: "hailgenetics/hail:0.2.126"
        memory: "256G"
        cpu: 64
        disks: "local-disk 1000 SSD"
    }

    output {
        File output_vcf = "~{output_path}"
    }
}

task plink2 {
    input {
        File vcf_file
        String output_prefix
        Int new_id_max_allele_len
    }

    command <<<
        set -e
        mkdir -p plink_output

        plink2 --vcf "~{vcf_file} \
        --make-pgen \
        --out plink_output/"~{output_prefix}" \
        --set-all-var-ids @:#\$r_\$a \
        --new-id-max-allele-len "~{new_id_max_allele_len} \
        --output-chr chrM \
        --chr 1-22
    >>>

    runtime {
        docker: "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
        memory: "16G"
        cpu: 4
        disks: "local-disk 100 SSD"
    }

    output {
        Directory plink_outputs = "plink_output"
    }
}
