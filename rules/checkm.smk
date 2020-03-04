checkm_data = config['checkm_data']

rule checkm_setup:
    output:
        join(checkm_data, 'selected_marker_sets.tsv')
    conda:
        "../envs/checkm.yaml"
    log:
       join(outdir, "logs/checkm/setup.log")
    shell:
        """

        wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

        tar -xzf checkm_data_2015_01_16.tar.gz -C {checkm_data}

        rm checkm_data_2015_01_16.tar.gz

        checkm data setRoot {checkm_data} 2> {log} 1>&2
        """


rule checkm_lineage_wf:
    input:
        setup=rules.checkm_setup.output,
        fasta=join(outdir, "assembled/spades/{sample}/contigs.fasta")
    output:
        lineage=join(outdir, "assembled/checkm/{sample}/output/lineage.ms"),
        cin=join(outdir, "assembled/checkm/{sample}/input/{sample}.fna")
    threads:
        8
    conda:
        "../envs/checkm.yaml"
    resources:
        mem_mb=12000
    log:
       join(outdir, "logs/checkm/{sample}.log")
    shell:
        """
        indir=$(dirname "{output.cin}")
        outdir=$(dirname "{output.lineage}")

        mkdir -p $indir
        cp {input.fasta} {output.cin}

        checkm lineage_wf --reduced_tree -t {threads} -x fna $indir $outdir 2> {log} 1>&2
        """


rule bwa_index:
    input:
        rules.checkm_lineage_wf.output.cin
    output:
        amb=join(outdir, "assembled/checkm/{sample}/input/{sample}.amb"),
        ann=join(outdir, "assembled/checkm/{sample}/input/{sample}.ann"),
        bwt=join(outdir, "assembled/checkm/{sample}/input/{sample}.bwt"),
        pac=join(outdir, "assembled/checkm/{sample}/input/{sample}.pac"),
        sa=join(outdir, "assembled/checkm/{sample}/input/{sample}.sa")
    log:
        join(outdir, "logs/bwa_index/{sample}.log")
    params:
        prefix=join(outdir, "assembled/checkm/{sample}/input/{sample}"),
        algorithm="bwtsw"
    wrapper:
        "0.35.0/bio/bwa/index"


rule bwa_mem:
    input:
        reads=[join(outdir, "trimmed/{sample}.trim.r1.fastq.gz"),
               join(outdir, "trimmed/{sample}.trim.r2.fastq.gz")],
        idx=rules.bwa_index.output.bwt
    output:
        join(outdir, "assembled/checkm/{sample}/input/{sample}.bam")
    log:
        join(outdir, "logs/bwa_mem/{sample}.log")
    params:
        index=join(outdir, "assembled/checkm/{sample}/input/{sample}"),
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.35.0/bio/bwa/mem"


rule samtools_index:
    input: 
        rules.bwa_mem.output
    output: 
        join(outdir, "assembled/checkm/{sample}/input/{sample}.bam.bai")
    params:
        "" # optional params string
    wrapper:
        "0.35.0/bio/samtools/index"


rule checkm_cov:
    input:
        binf=rules.checkm_lineage_wf.output.cin,
        bam=rules.bwa_mem.output,
        bai=rules.samtools_index.output
    output:
        join(outdir, "assembled/checkm/{sample}/output/coverage.tsv")
    conda:
        "../envs/checkm.yaml"
    log:
        join(outdir, "logs/checkm/{sample}.cov.log")
    shell:
        """
        bindir=$(dirname "{input.binf}")
        checkm coverage $bindir {output} {input.bam}
        """


rule checkm_tetra:
    input:
        rules.checkm_lineage_wf.output.cin
    output:
        join(outdir, "assembled/checkm/{sample}/output/tetra.tsv")
    conda:
        "../envs/checkm.yaml"
    log:
        join(outdir, "logs/checkm/{sample}.tetra.log")
    shell:
        """
        checkm tetra {input} {output}
        """


rule checkm_dist_plot:
    input:
        tetra=rules.checkm_tetra.output,
        outf=rules.checkm_lineage_wf.output.lineage,
        binf=rules.checkm_lineage_wf.output.cin
    output:
        join(outdir, "assembled/checkm/{sample}/output/plots/{sample}.ref_dist_plots.png")
    conda:
        "../envs/checkm.yaml"
    log:
        join(outdir, "logs/checkm/{sample}.dist_plots.log")
    shell:
        """
        bindir=$(dirname "{input.binf}")
        plotdir=$(dirname "{output}")
        outdir=$(dirname "{input.outf}")

        checkm dist_plot $outdir $bindir $plotdir {input.tetra} 95
        """


rule checkm_map:
    input:
        expand(join(outdir, "assembled/checkm/{sample}/input/{sample}.bam"), sample=samples.index)


rule checkm_all:
    input:
        expand(join(outdir, "assembled/checkm/{sample}/output/lineage.ms"), sample=samples.index),
        expand(join(outdir, "assembled/checkm/{sample}/output/plots/{sample}.ref_dist_plots.png"), sample=samples.index),
