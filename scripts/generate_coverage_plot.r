library("biomaRt")
library(tidyverse)
library(argparse)
library(stringr)
library(yaml)



# Argument parsing
parser = ArgumentParser()
parser$add_argument("coverage_tables", help = paste("Comma separated string containing paths to coverage tables. ",
                                                    "tables are named of form <gene_id>_<condition>.bedgraph",
                                                  sep = ""))
parser$add_argument("outdir", help="Path to output directory")
parser$add_argument("config", help="Path to config file")
parser$add_argument("script_dir", help="path to directory containing scripts")
parser$add_argument("transcripts")
parser$add_argument("symbol")
# parser$add_argument("-l", "--location", help = "Location of the pipeline folder", default = "D:/work/rna-seq_pipelines/scripts/analysis")
args = parser$parse_args()

# Testing parameters
# {
#   args = list(coverage_tables="output/gene_genomecovs/ENSMUSG00000030795_kd-aso_coverage.bedgraph,output/gene_genomecovs/ENSMUSG00000030795_control_coverage.bedgraph,output/gene_genomecovs/ENSMUSG00000030795_kd-control_coverage.bedgraph,output/gene_genomecovs/ENSMUSG00000030795_homo_coverage.bedgraph,output/gene_genomecovs/ENSMUSG00000030795_hetero_coverage.bedgraph",
#               outdir="output/coverage_figures",
#               config="D:/work/sc/workflows/runs/fus_coverage/config.yaml",
#               script_dir="D:/work/sc/workflows/generate_coverage/scripts",
#               transcripts = "ENSMUST00000123151,ENSMUST00000077609,ENSMUST00000121616,ENSMUST00000174196,ENSMUST00000141997,ENSMUST00000106251",
#               symbol="Fus")
# }
source(file.path(args$script_dir, "coverage.r"))



  symbol = args$symbol
  tables = str_split(args$coverage_tables, ",")[[1]]
  data = data.frame(str_match(tables, "([\\w-]+)_([\\w-]+)_coverage"))
  colnames(data) = c("name", "gene", "condition")
  gene = data$gene[[1]] %>% as.character()
  data$paths = tables
  print(args$transcripts)
  transcripts = args$transcripts %>% strsplit(",") %>% unlist()
  config= read_yaml(file=args$config)
  # order of the conditions for plotting is assumed to be the order of cond_map values
  condition_order = config$cond_map %>% unlist() %>% unique()
  condition_colors = config$cond_colors
# Read and relabel the genomecov data
genomecovs = lapply(data$paths, function(x){
  a = read_tsv(x, col_names = F)
  colnames(a) = c("chromosome", "start", "end", paste("s", 1:(ncol(a) - 3), sep=""))
  a = norm_chrom(a)
  a
  })
names(genomecovs) = data$condition
if (any(names(genomecovs) %>% is.na)){
  print(genomecovs)
  stop("Name with NA detected")
}
# print(data)
# print(genomecovs)
# stop()


print("CALCULATING LIMITS")
ylim = c(0, max(lapply(genomecovs, get_limits) %>% unlist))



print("getting ensembl object")
# Get the biomart ensembl object for coverage generation
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
print("Generating biomTrack")
print(paste("Variables", symbol, transcripts))
# Generate the Biomart track
# biomTrack <- BiomartGeneRegionTrack(genome="GRCm38.p5",
#                                     biomart=ensembl, stacking = "squish",
#                                     filters = list(ensembl_transcript_id=c(transcripts)))
# stop()
biomTrack <- BiomartGeneRegionTrack(genome="GRCm38.p5",
                                    # symbol=toupper(symbol), # Can't do b oth symbol and transcript filter
                                    name=paste(symbol, "\nstructure", sep=""),
                                    biomart=ensembl, stacking = "squish",
                                    background.title = "transparent",
                                    col = "black",
                                    fill = "black",
                                    col.line = "black",
                                    protein_coding = "black",
                                    utr3 = "black",
                                    utr5 = "black",
                                    non_coding = "black",
                                    min.height = 0,
                                    filters = list(ensembl_transcript_id=c(transcripts))
                                    )

# Generate the genomecov tracks
print("GENERATING GENOEMCOV TRACKS")
generate_datatrack_wrapper = function(index, ylim) {
  # Wrapper around datatrack generation to make customization of appearance easier
  
  genomecov_data = genomecovs[[index]]
  genomecov_name = names(genomecovs)[[index]]
  print(index)
  print(genomecov_name)
  col = ifelse(!is.na(condition_colors[[genomecov_name]]),
               condition_colors[[genomecov_name]],
               "black")
  print(col)
  symbol = "Fus"
  make_datatrack(genomecov_data, name=genomecov_name, col=col, symbol=symbol, ylim=ylim)
}
datatracks = lapply(seq_along(genomecovs), generate_datatrack_wrapper, ylim=ylim)
names(datatracks) = names(genomecovs)
print(datatracks)
print("GENERATING IMAGES")

png(filename = file.path(args$outdir, paste(gene, "_coverage.png", sep="")), width = 4, height=6, units = "in", res=300)
plotTracks(c(biomTrack, datatracks[condition_order]))
dev.off()

svg(filename = file.path(args$outdir, paste(gene, "_coverage.svg", sep="")), width = 4, height=6)
plotTracks(c(biomTrack, datatracks[condition_order]))
dev.off()
