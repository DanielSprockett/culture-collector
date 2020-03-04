gtdb_data = config['gtdb_data']

rule gtdb_setup:
    output:
        join(gtdb_data, '89.0', 'bac120_taxonomy_r89.tsv')
    conda:
        "../envs/gtdb.yaml"
    log:
       join(outdir, "logs/gtdb/setup.log")
    shell:
        """

        wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
        tar -xvzf gtdbtk_r89_data.tar.gz -C {gtdb_data}

        rm xvzf gtdbtk_r89_data.tar.gz
        """

rule taxonomy_gtdb:
    input:
        fasta=rules.assemble_spades.output
    output:
        touch(join(outdir, "taxonomy/gtdb/{sample}/gtdbtk.done"))
    log:
        join(outdir, "logs/gtdb/{sample}.log")
    conda:
        "../envs/gtdb.yaml"
    params:
        config["params"]["gtdb"]
    shell:
        """
        export GTDBTK_DATA_PATH={gtdb_data}

        outdir=$(dirname "{output[0]}")

        echo -e "{input[0]}\t{wildcards.sample}" > $outdir/file.txt

        gtdbtk classify_wf --batchfile $outdir/file.txt --out_dir $outdir 2> {log} 1>&2
        """

