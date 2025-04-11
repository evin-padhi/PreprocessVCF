version 1.0

workflow VDS_to_VCF_pipeline {
    input {
        String vds_url
        String rnaseq_samples_tsv
        String identifier
        String output_bucket
    }

    call ConvertVdsToDenseMt {
        input:
            vds_url = vds_url,
            rnaseq_samples_tsv = rnaseq_samples_tsv,
            identifier = identifier,
            output_bucket = output_bucket,
    }

    call SplitMultiAllelic {
        input:
            input_mt_path = ConvertVdsToDenseMt.output_mt_path,
            identifier = identifier,
            output_bucket = output_bucket
    }

    call ExportVcf {
        input:
            input_mt_path = SplitMultiAllelic.split_mt_path,
            identifier = identifier,
            output_bucket = output_bucket
    }

    output {
        File final_vcf_bgz = ExportVcf.output_vcf_bgz
    }
}

task ConvertVdsToDenseMt {
    input {
        String vds_url
        String rnaseq_samples_tsv
        String output_bucket
        String identifier

        String docker = "quay.io/jonnguye/hail:latest"
        Int memory = 512
        Int cpu = 64
        Int disk_size = 500
        string disk_type = SSD
    }

    command {
        curl -O https://raw.githubusercontent.com/jonnguye/PreprocessVCF/NotebookToWDL/ConvertVDSToDenseMT.py
        python3 ConvertVDSToDenseMT.py ~{vds_url} ~{rnaseq_samples_tsv} ~{output_bucket} ~{identifier}
    }

    runtime {
        docker: ~{docker}
        memory: "~{memory} GB"
        cpu: ~{cpu}
        disks: "local-disk ~{disk_size} ~{disk_type}"
    }

    output {
        String output_mt_path = "~{output_bucket}/~{identifier}.mt"
    }
}

task SplitMultiAllelic {
    input {
        String input_mt_path
        String output_bucket
        String identifier

        String docker = "quay.io/jonnguye/hail:latest"
        Int memory = 512
        Int cpu = 64
        Int disk_size = 500
        string disk_type = SSD
    }

    command {
        curl -O https://raw.githubusercontent.com/jonnguye/PreprocessVCF/NotebookToWDL/SplitMultiAllelic.py
        python3 SplitMultiAllelic.py ~{input_mt_path} ~{output_bucket} ~{identifier}
    }

    runtime {
        docker: ~{docker}
        memory: "~{memory} GB"
        cpu: ~{cpu}
        disks: "local-disk ~{disk_size} ~{disk_type}"
    }

    output {
        String split_mt_path = "~{output_bucket}/~{identifier}_split.mt"
    }
}

task ExportVcf {
    input {
        String input_mt_path
        String output_bucket
        String identifier
        
        String docker = "quay.io/jonnguye/hail:latest"
        Int memory = 256
        Int cpu = 32
        Int disk_size = 500
        string disk_type = SSD
    }

    command {
        curl -O https://raw.githubusercontent.com/jonnguye/PreprocessVCF/NotebookToWDL/ExportVCF.py
        python3 ExportVCF.py ~{input_mt_path} ~{output_bucket} ~{identifier}
    }

    runtime {
        docker: ~{docker}
        memory: "~{memory} GB"
        cpu: ~{cpu}
        disks: "local-disk ~{disk_size} ~{disk_type}"
    }

    output {
        File output_vcf_bgz = "~{output_bucket}/~{identifier}.vcf.bgz"
    }
}
