generate_coverage
Generates coverage plots for a given set of genes and conditions from a set of bamfiles
Input required:
	bamfiles in the form of <sample_id>

Data preparation:
	Configuration of the run:
		Create a config.yaml file in the run directory (working directory) with the following parameters
			cond_map: mapping of sample names to condition names
			indir: Input directory containing bamfiles to be used. Each bamfile is expected to have name <sample_id>.bam
			target_genes: Gene for which to generate coverage for
			gene_coords: Path to gene_coord annotations



Notes:
	Bedtools needed, thus run from within a bash shell



