# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
import pandas as pd
from snakemake.utils import validate, min_version
from os.path import join

#### set minimum Snakemake version
min_version("5.9.1")

#### load config and sample sheets

configfile: "config.yaml"
validate(config, schema='schemas/config.schema.yaml')

outdir = config['outdir']

samples = pd.read_csv(config['samples'],
                      sep='\t',
                      comment='#',
                      dtype={'sample': str}).set_index('sample', drop=False)
validate(samples, schema='schemas/samples.schema.yaml')

units = pd.read_csv(config['units'], sep='\t', comment='#').set_index('sample', drop=False)
# units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema='schemas/units.schema.yaml')

report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

# Specify location for temp files, e.g. scrach space:
TMP_DIR_ROOT = config['temp_dir']


rule all:
    input:
        join(outdir, "qc/multiqc.html"),
        join(outdir, "assembled/multiqc.html"),
        expand(join(outdir, "assembled/checkm/{sample}/output/lineage.ms"), sample=samples.index),
        # expand(join(outdir, "assembled/checkm/{sample}/output/plots/cov_pca.done"), sample=samples.index),
        expand(join(outdir, "assembled/checkm/{sample}/output/plots/{sample}.ref_dist_plots.png"), sample=samples.index),
        # expand(join(outdir, "assembled/checkm/{sample}/output/plots/{sample}.paralel_coord_plot.png"), sample=samples.index),
        # expand(join(outdir, "assembled/checkm/{sample}/output/plots/bin_qa_plot.png"), sample=samples.index),
        expand(join(outdir, "taxonomy/gtdb/{sample}/{sample}.bac120.summary.tsv"), sample=samples.index)

        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.


include: "rules/qc.smk"
include: "rules/assemble.smk"
include: "rules/checkm.smk"
include: "rules/taxonomy.smk"