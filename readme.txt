generate_coverage

Data preparation:
	Configuration of the run:
		Create a config.yaml file in the run directory (working directory) with the following parameters
			cond_map: mapping of sample names to condition names
			indir: Input directory containing bamfiles to be used. Each bamfile is expected to have name <sample_id>.bam
			target_genes: Gene for which to generate coverage for
			gene_coords: Path to gene_coord annotations



Generates coverage plots given:
	A set of bams
	Condition annotations for each bam (Table containing sample name, file_name, condition)
	A set of genes and the selected transcripts

Notes:
	Bedtools needed, thus run from within a bash shell



