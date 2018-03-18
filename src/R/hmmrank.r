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

SCRIPT_VERSION = '1.0.1'

# Get arguments
option_list = list(
  make_option(
    c('--minscore'), type='double', default=0.0,
    help='Minimum score, default %default'
  ),
  make_option(
    c('--outfile'), type='character',
    help='Name of output file'
  ),
  make_option(
    c('--scorefile'), type='character',
    help='Name of tsv file containing accept scores (GA in score_type): profile, score_type and seq_score'
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

logmsg(sprintf("hmmrank.r %s", SCRIPT_VERSION))

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
  if ( t %>% nrow() > 0 ) tblout <- union(tblout, t %>% select(accno, profile, evalue, score))
}

# If we're given a scorefile as an option, read and left join
if ( length(opt$options$scorefile) > 0 ) {
  logmsg(sprintf("Joining in scores from %s, setting minimum score for absents to %5.2f", opt$options$scorefile, opt$options$minscore))
  tblout <- tblout %>%
    left_join(
      read_tsv(
        opt$options$scorefile,
        col_types = cols(.default = col_character(), seq_score = col_double(), domain_score = col_double())
      ) %>%
        filter(score_type == 'GA'),
      by = 'profile'
    ) %>%
    replace_na(list('seq_score' = opt$options$minscore))
} else {
  logmsg(sprintf("Setting minimum scores to %5.2f", opt$options$minscore))
  tblout <- tblout %>% mutate(seq_score = opt$options$minscore)
}

# Calculate rank
logmsg("Calculating ranks")
tblout <- tblout %>% 
  filter(score >= seq_score) %>%
  group_by(accno) %>% mutate(rank = rank(desc(score)))

logmsg(sprintf("Writing to %s", opt$options$outfile))
write_tsv(
  tblout %>% select(accno, profile, rank, evalue, score) %>% arrange(accno, rank),
  opt$options$outfile
)

logmsg("Done")
