version 1.0

workflow FilterAndWriteMTWorkflow {
    input {
        String matrix_table
        String samples_table
        String output_checkpoint
    }

    call FilterAndWriteMTTask {
        input:
            matrix_table = matrix_table,
            samples_table = samples_table,
            output_checkpoint = output_checkpoint
    }

    output {
        String checkpoint_path = FilterAndWriteMTTask.output_checkpoint
    }
}

task FilterAndWriteMTTask {
    input {
        String matrix_table
        String samples_table
        String output_checkpoint
    }

    command <<<
        set -e

        # Download the Python script
        curl -O https://raw.githubusercontent.com/jonnguye/PreprocessVCF/NotebookToWDL/filter_and_write_mt.py

        # Run the Python script
        python3 filter_and_write_mt.py \
            --matrix_table "~{matrix_table}" \
            --samples_table "~{samples_table}" \
            --output_checkpoint "~{output_checkpoint}"
    >>>

    runtime {
        docker: "quay.io/jonnguye/hail:latest"
        memory: "512G"
        cpu: 64
        disks: "local-disk 1000 SSD"
    }

    output {
        String output_checkpoint = "~{output_checkpoint}"
    }
}
