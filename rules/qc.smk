def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()


rule file_rename:
    input:
        get_fastq
    output:
        r1=temp(join(outdir, "input/{sample}.r1.fastq.gz")),
        r2=temp(join(outdir, "input/{sample}.r2.fastq.gz"))
    run:
        shell("""
              cp {input[0]} {output.r1}
              cp {input[1]} {output.r2}
              """)


rule cutadapt_pe:
    input:
        rules.file_rename.output
    output:
        fastq1=join(outdir, "trimmed/{sample}.trim.r1.fastq.gz"),
        fastq2=join(outdir, "trimmed/{sample}.trim.r2.fastq.gz"),
        qc=join(outdir, "trimmed/{sample}.qc.txt")
    params:
        config["params"]["cutadapt-pe"]
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters = config["params"]["cutadapt-pe"]["adapters"],
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        others = config["params"]["cutadapt-pe"]["others"]
    log:
        join(outdir, "logs/cutadapt/{sample}.log")
    wrapper:
        "0.49.0/bio/cutadapt/pe"


rule fastqc_r1:
    input:
        rules.file_rename.output.r1
    output:
        html=join(outdir, "qc/fastqc/{sample}.r1_fastqc.html"),
        zip=join(outdir, "qc/fastqc/{sample}.r1_fastqc.zip")
    params: ""
    log:
        join(outdir, "logs/fastqc/{sample}.r1.log")
    wrapper:
        "0.49.0/bio/fastqc"


rule fastqc_r2:
    input:
        rules.file_rename.output.r2
    output:
        html=join(outdir, "qc/fastqc/{sample}.r2_fastqc.html"),
        zip=join(outdir, "qc/fastqc/{sample}.r2_fastqc.zip")
    params: ""
    log:
        join(outdir, "logs/fastqc/{sample}.r2.log")
    wrapper:
        "0.49.0/bio/fastqc"


rule fastqc_trim_r1:
    input:
        rules.cutadapt_pe.output['fastq1']
    output:
        html=join(outdir, "trimmed/fastqc/{sample}.r1_fastqc.html"),
        zip=join(outdir, "trimmed/fastqc/{sample}.r1_fastqc.zip")
    params: ""
    log:
        join(outdir, "logs/fastqc/{sample}.trim.r1.log")
    wrapper:
        "0.49.0/bio/fastqc"


rule fastqc_trim_r2:
    input:
        rules.cutadapt_pe.output['fastq2']
    output:
        html=join(outdir, "trimmed/fastqc/{sample}.r2_fastqc.html"),
        zip=join(outdir, "trimmed/fastqc/{sample}.r2_fastqc.zip")
    params: ""
    log:
        join(outdir, "logs/fastqc/{sample}.trim.r2.log")
    wrapper:
        "0.49.0/bio/fastqc"


rule multiqc:
    input:
        expand(rules.fastqc_r1.output, sample=samples.index, unit=units.index),
        expand(rules.fastqc_trim_r1.output, sample=samples.index, unit=units.index),
        expand(rules.fastqc_r2.output, sample=samples.index, unit=units.index),
        expand(rules.fastqc_trim_r2.output, sample=samples.index, unit=units.index),
        expand(rules.cutadapt_pe.output.qc, sample=samples.index, unit=units.index)
    output:
        join(outdir, "qc/multiqc.html")
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        join(outdir, "logs/multiqc.log")
    wrapper:
        "0.49.0/bio/multiqc"