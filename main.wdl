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

    output {
        File output_vcf = WriteVCFTask.output_vcf
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
