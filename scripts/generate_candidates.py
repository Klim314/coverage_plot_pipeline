"""
Subsets a gene annotation table containing gene/transcrpt start/ends
"""

def create_subset()


if __name__ == "__main__":
    gene_data = "/mnt/d/work/sc/data/biomart_dumps/ensembl_gene_start_end.tsv"
    gene_data = pd.read_csv(gene_data, sep="\t")
    gene_data.columns = ["ens_gene", "ens_trans", "chrom", 'start', 'end', 'tsl', 'symbol']
    gene_data = gene_data[gene_data["ens_gene"].isin(target_genes)]
    gene_data.to_csv(pjoin(output_dir, "gene_table_all-tsl.tsv"), sep="\t", index=None)
    gene_data = gene_data[(gene_data['tsl'].str.startswith("tsl1")) | (gene_data['tsl'].str.startswith("tsl2"))]
    gene_data.to_csv(pjoin(output_dir, "gene_table.tsv"), sep="\t", index=None)