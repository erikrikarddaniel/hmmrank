#!/usr/bin/env Rscript

# hmmrank.r
#
# Ranks sequences according to their score against different hmm profiles.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Get arguments
option_list = list(
  make_option(
    c('--outfile'), type='character',
    help='Name of output file'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] file0.tblout .. filen.tblout", 
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

# Table to fill with data from files
tblout <- tibble(
  accno = character(), profile = character(),
  evalue = double(), score = double()
)

# Read all the tblout files
for ( tbloutfile in grep('\\.tblout', opt$args, value=TRUE) ) {
  logmsg(sprintf("Reading %s", tbloutfile))
  t <-  read_fwf(
    tbloutfile, fwf_cols(content = c(1, NA)), 
    col_types = cols(content = col_character()), 
    comment='#'
  ) %>% 
    separate(
      content, 
      c('accno', 't0', 'profile', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'), 
      '\\s+', 
      extra='merge',
      convert = T
    )
  tblout <- union(tblout, t %>% select(accno, profile, evalue, score))
}

# Calculate rank
tblout <- tblout %>% group_by(accno) %>% mutate(rank = rank(desc(score)))

write_tsv(
  tblout %>% select(accno, profile, rank, evalue, score) %>% arrange(accno, rank),
  opt$options$outfile
)

logmsg("Done")
