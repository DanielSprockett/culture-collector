rule assemble_spades:
    input:
        fastq1=join(outdir, "trimmed/{sample}.trim.r1.fastq.gz"),
        fastq2=join(outdir, "trimmed/{sample}.trim.r2.fastq.gz")
    output:
        join(outdir, "assembled/spades/{sample}/contigs.fasta")
    log:
        join(outdir, "logs/spades/{sample}.log")
    conda:
        "../envs/assemble.yaml"
    params:
        config["params"]["spades"]
    resources:
        mem_mb=9000
    threads:
        4
    shell:
        """
        m_gb=$(( {resources.mem_mb} / 1000 ))

        outdir=$(dirname "{output[0]}")

        spades.py -t {threads} -m $m_gb \
        -1 {input.fastq1} -2 {input.fastq2} -o $outdir {params} \
        2> {log} 1>&2
        """


rule assemble_quast:
    input:
        rules.assemble_spades.output
    output:
        join(outdir, "assembled/quast/{sample}/report.tsv")
    log:
        join(outdir, "logs/quast/{sample}.log")
    conda:
        "../envs/assemble.yaml"
    params:
        config["params"]["quast"]
    shell:
        """
        outdir=$(dirname "{output[0]}")

        quast.py -t {threads} -o $outdir/../ \
        {input} 2> {log} 1>&2
        """



rule assemble_metaqc:
    input:
        expand(rules.assemble_quast.output, sample=samples.index, unit=units.index)
    output:
        join(outdir, "assembled/multiqc.html")
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        join(outdir, "logs/assemble_multiqc.log")
    wrapper:
        "0.49.0/bio/multiqc"