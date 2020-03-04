rule taxonomy_gtdb:
    input:
        rules.assemble_spades.output
    output:
        join(outdir, "taxonomy/gtdb/{sample}/{sample}.bac120.summary.tsv")
    log:
        join(outdir, "logs/gtdb/{sample}.log")
    conda:
        "../envs/gtdb.yaml"
    params:
        config["params"]["gtdb"]
    shell:
        """
        outdir=$(dirname "{output[0]}")

        echo -e "{input[0]}\t{wildcards.sample}" > $outdir/file.txt

        gtdbtk classify_wf --batchfile $outdir/file.txt --out_dir $outdir 2> {log} 1>&2
        """

