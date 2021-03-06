from collections import defaultdict
import os
import pandas as pd

# Load the snakemake configfile
configfile: "config.yaml"


def pjoin(*args):
    return "/".join(args)

# Constants and derived variables

input_dir = "input"
output_dir = "output"
temp_dir = "temp"
for i in [input_dir, output_dir, temp_dir]:
    if not os.path.exists(i):
        os.makedirs(i)

# Process the condition map
cond_map = config["cond_map"]
conditions = {item for key, item in cond_map.items()}  # Set of all conditions in the experiment
cond_sample = defaultdict(list)
for sample, condition in cond_map.items():
    cond_sample[condition].append(sample)
    cond_sample[condition].sort()



target_genes = config["target_genes"]

# Take in annotations for the gene coordinate data and normalize their names
# End result is a table containing only Transcript Support Leel (1) data and appris primary data
gene_data = config["gene_coords"]
gene_data = pd.read_csv(gene_data, sep="\t")
gene_data = gene_data.rename(columns={"Gene stable ID": "ens_gene",
                                      "Transcript stable ID": "ens_trans",
                                      "Chromosome/scaffold name": "chrom",
                                      "Gene start (bp)": "start",
                                      "Gene end (bp)": "end",
                                      "Transcript support level (TSL)": "tsl",
                                      "APPRIS annotation": "appris",
                                      "Gene name": "symbol"})

gene_data = gene_data[gene_data["ens_gene"].isin(target_genes)]
gene_data.to_csv(pjoin("output", "gene_table_all-tsl.tsv"), sep="\t", index=None)
gene_data = gene_data[gene_data['tsl'].str.startswith("tsl1") |
                      gene_data['tsl'].str.startswith("tsl2") |
                      gene_data['appris'].str.startswith("principal")]
gene_data.to_csv(pjoin("output", "gene_table.tsv"), sep="\t", index=None)

# Structure the gene range data in a nicer form
g_ranges = dict()

for row in gene_data.iterrows():
    index, data = row
    ens_gene, ens_trans, chrom, start, end, tsl, appris, symbol  = data.tolist()
    if ens_gene in target_genes:
        g_ranges[ens_gene] = (chrom, start, end)




# This is the expected outputs, a coverage bedgraph for each gene and condition
# print(expand(pjoin("output", "gene_genomecovs", "{gene}_{cond}_coverage.bedgraph"),
#                cond=conditions,
#                gene=target_genes))
print(expand(pjoin("output", "coverage_figures", "{gene_id}_coverage.png"), gene_id=target_genes))

rule all:
    input:
        expand(pjoin("output", "coverage_figures", "{gene_id}_coverage.png"), gene_id=target_genes)
        # expand(pjoin("output", "gene_genomecovs", "{gene}_{cond}_coverage.bedgraph"),
        #        cond=conditions,
        #        gene=target_genes)
        # expand(pjoin("output", "summarized_genomecovs", "{cond}_summarized.bedgraph"), cond=conditions)

rule prepare_gene_beds:
    """Generates bedifles that are essentially just the start and end of the bed
    """
    output: "temp/{gene}.bed"
    run: 
        ens_id = wildcards.gene
        with open(str(output), "w") as f:
            f.write("{}\t{}\t{}\n".format(*g_ranges[ens_id]))

rule genomecov:
    input: "/".join([input_dir, "{sample}.bam"])
    output: pjoin("output", "genomecovs", "{sample}_cov.bedgraph")
    run:
        # print(output)
        out_dir = os.path.split(str(output))[0]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        shell("bedtools genomecov -ibam {input} -bg -split > {output}")

rule summarize_genomecovs:
    input:
        bedgraphs = lambda wildcards: expand(pjoin("output", "genomecovs", "{sample}_cov.bedgraph"),
                                             sample=cond_sample[wildcards.cond])
    output: pjoin("output", "summarized_genomecovs", "{cond}_summarized.bedgraph")
    params:
        outdir = pjoin("output", "summarized_genomecovs")
    run:
        if not os.path.exists(params.outdir):
            os.makedirs(params.outdir)
        # UnionBedG will fail if only one
        if len(input.bedgraphs) == 1:
            bedgraph = input.bedgraphs[0]
            shell("cp {bedgraph} {output}")
        else:
            shell("bedtools unionbedg -i {input.bedgraphs} > {output}")



rule generate_subsets:
    """
    generates the various subbedgraphs that are desired
    specifically, we're grabbing the bedgraphs for specific gene regions
    """
    input: 
        bedgraph = pjoin("output", "summarized_genomecovs", "{cond}_summarized.bedgraph"),
        gene_bed = "temp/{gene}.bed"
    output: pjoin("output", "gene_genomecovs", "{gene}_{cond}_coverage.bedgraph")

    run:
        out_dir = os.path.split(str(output))[0]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        print("Generating gene: {}, condition: {}".format(wildcards.gene, wildcards.cond))
        shell("bedtools intersect -a {input.bedgraph} -b {input.gene_bed} > output/gene_genomecovs/{wildcards.gene}_{wildcards.cond}_coverage.bedgraph")


print(expand(pjoin("output", "gene_genomecovs", "{{gene_id}}_{cond}_coverage.bedgraph"), cond=conditions))
rule generate_figures:
    """ Generates genomecoverage figures compiling all conditions into a single plot
    """
    input: expand(pjoin("output", "gene_genomecovs", "{{gene_id}}_{cond}_coverage.bedgraph"), cond=conditions)
    output:
        png = pjoin("output", "coverage_figures", "{gene_id}_coverage.png"),
        svg = pjoin("output", "coverage_figures", "{gene_id}_coverage.svg"),
    params:
        outdir = pjoin("output", "coverage_figures"),
        script_dir = srcdir("scripts"),
        coverage_script = srcdir("scripts/generate_coverage_plot.r")
    run:
        if not os.path.exists(params.outdir):
            os.makedirs(params.outdir)

        table_paths = ",".join(input)
        # Indentify transcripts of interest for the plot and create a nice little string to feed in later
        reduced_gene_data = gene_data.loc[gene_data["ens_gene"] == wildcards.gene_id]
        print(reduced_gene_data)
        transcripts = list(set(reduced_gene_data["ens_trans"]))
        transcripts = ",".join(transcripts)
        symbol = reduced_gene_data["symbol"].tolist()[0]

        shell("Rscript {params.coverage_script} {table_paths} {params.outdir} config.yaml {params.script_dir} {transcripts} {symbol}")

# rule generate_figures:
#     """
#     Generates the figures using a simple R script
#     """
#     input:
#         gene_table = pjoin(temp_dir, "gene_table.tsv"),
#         bedgraphs = expand(pjoin("output", 'gene_genomecovs', '{gene}_dko_coverage.bedgraph'), target_genes)
#     output: expand(pjoin("output", "figures", "{gene}_figure.png"), gene=target_genes)
#     run:

