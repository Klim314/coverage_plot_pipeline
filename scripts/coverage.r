library("data.table")
library("biomaRt")
library("Gviz")
library("GenomicRanges")
library("GenomicAlignments")
library("dplyr")

norm_chrom = function(table){
  # Normalizes chromosome names to UCSC format
  table = table
  table = table %>%
    mutate(chromosome = ifelse(chromosome == 'X', 'chrX', chromosome)) %>% 
    mutate(chromosome = ifelse(chromosome == 'Y', 'chrY', chromosome))
  table = table %>% filter(chromosome %in% c(1:200, 'chrX', 'chrY'))
  return(table)
}

make_datatrack = function(table, col, name, symbol, ylim=c(0, 150), col.axis = "black"){
  # Creates a histogram plot for a given table of genomecoverage data
  # args
  #   table: table object containing the following columns
  #   col: Color of 
  #   name: Name to use
  #   symbol: Symbol of gene to investigate
  #   ylm: Y-axis limits for plotting
  DataTrack(range = table,
            genome="GRCm38.p5", 
            symbol = symbol,
            name=name,
            biomart = ensembl,
            type = 'hist',
            ylim = ylim, 
            #            window = window_size,
            col.histogram = col,
            fill.histogram = col,
            background.title = "transparent",
            col.id="#000000",
            col.title = "black",
            col.axis = col.axis,
            rot.title = 0
  )
}

make_cliptrack = function(table, col, name, symbol, col.axis = "black"){
  # Constructs a datatrack containing CLIP peaks represented as single lines
  #
  #
  #
  DataTrack(range = table,
            genome="GRCm38.p5", 
            symbol = symbol,
            name=name,
            biomart = ensembl,
            type = 'hist',
            windowSize = 1,
            #            window = window_size,
            col.histogram = col,
            fill.histogram = col, 
            background.title = "transparent",
            col.id="#000000",
            col.title = "black",
            col.axis = col.axis,
            rot.title = 0
            
  )
}

get_limits = function(table){
  # Determines the Y-axis limits from a genomecoverage bedfile
  # First three columns must be Chr, start, end
  # Takes the mean of all remaining columns
  return(
    c(0,
      table[, 4:ncol(table)] %>%
        rowMeans %>%
        max %>%
        as.integer
    )
  )
}
