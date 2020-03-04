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

rule gtdb_prep:
    input:
        fastas=expand(rules.assemble_quast.output, sample=samples.index, unit=units.index)
    output:
        join(outdir, "taxonomy", "gtdb", "batchfile.txt")
    log:
        join(outdir, "logs", "gtdb", "prep.log")
    run:
        with open(output[0], 'w') as f:
            for s in samples.index:
                fasta = join(outdir, "assembled/spades/%s/contigs.fasta" % s)
                if fasta not in input.fastas:
                    raise ValueError("Can't find input file %s" % fasta)
                f.write("{0}\t{1}".format(fasta, s))


rule taxonomy_gtdb:
    input:
        rules.gtdb_prep.output
    output:
        touch(join(outdir, "taxonomy/gtdb/gtdbtk.done"))
    log:
        join(outdir, "logs/gtdb/batch.log")
    conda:
        "../envs/gtdb.yaml"
    params:
        config["params"]["gtdb"]
    resources:
        mem_mb=200000
    threads:
        32
    shell:
        """
        export GTDBTK_DATA_PATH={gtdb_data}

        outdir=$(dirname "{output[0]}")

        gtdbtk classify_wf --cpus {threads} {params} --batchfile $outdir/file.txt --out_dir $outdir 2> {log} 1>&2
        """

