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
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

SCRIPT_VERSION = '1.3.0'

# Get arguments
# For interactive testing:
# opt <- list('options' = list('verbose' = TRUE, qfromfname = FALSE, minscore = 0, scorefile = 'hmmrank.02.profile_scores.tsv'), 'args' = c('hmmrank.02.d/NrdAe.tblout', 'hmmrank.02.d/NrdAg.tblout', 'hmmrank.02.d/NrdAh.tblout', 'hmmrank.02.d/NrdAi.tblout', 'hmmrank.02.d/NrdAk.tblout', 'hmmrank.02.d/NrdAm.tblout', 'hmmrank.02.d/NrdAn.tblout', 'hmmrank.02.d/NrdAq.tblout', 'hmmrank.02.d/NrdA.tblout', 'hmmrank.02.d/NrdAz3.tblout', 'hmmrank.02.d/NrdAz4.tblout', 'hmmrank.02.d/NrdAz.tblout'))
#
# Annotation table testing:
# opt <- list('options' = list('verbose' = TRUE, qfromfname = TRUE, minscore = 0, scorefile = 'hmmrank.06.profile_scores.tsv', annottable = 'hmmrank.06.annottable.tsv'), 'args' = c('hmmrank.06.d/arCOG01044.tblout', 'hmmrank.06.d/ENOG4108Z73.tblout', 'hmmrank.06.d/PF00590.tblout', 'hmmrank.06.d/PF06180.tblout', 'hmmrank.06.d/TIGR01466.tblout', 'hmmrank.06.d/TIGR01467.tblout', 'hmmrank.06.d/TIGR01469.tblout'))
option_list = list(
  make_option(
    '--annottable', type = "character",
    help = "Name of table with annotation assignments, at a minium must contain 'protein' and 'profile' columns. The 'protein' column needs to be unique. The 'profile' column might contain multiple profiles, separated by '&'."
  ),
  make_option(
    c('--minscore'), type='double', default=0.0,
    help='Minimum score, default %default'
  ),
  make_option(
    c("--only_best_scoring"), action="store_true", default=FALSE, 
    help="Output only the best scoring, i.e. scores < 2"
  ),
  make_option(
    c('--outfile'), type='character',
    help='Name of output file'
  ),
  make_option(
    "--qfromfname", action = "store_true", default = FALSE, 
    help = "Extract query name from file name -- basename until first period -- instead of taking it from the query column in each data file, default false."
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

# Read all the tblout files
tlist <- list()
i <- 1
for ( tbloutfile in grep('\\.tblout', opt$args, value=TRUE) ) {
  p <- basename(tbloutfile) %>% str_remove('\\..*')
  logmsg(sprintf("Reading %s", tbloutfile))
  t <- read_fwf(
    tbloutfile, fwf_cols(content = c(1, NA)), 
    col_types = cols(content = col_character()), 
    comment='#'
  ) %>% 
    separate(
      content, 
      c('accno', 't0', 'profile', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'), 
      '\\s+',  extra='merge', convert = FALSE
    ) %>%
    transmute(accno, profile, evalue = as.double(evalue), score = as.double(score)) %>%
    data.table()
  if ( opt$options$qfromfname ) t$profile <- p
  if ( nrow(t) > 0 ) {
    tlist[[i]] <- t
    i <- i + 1
  }
}
if ( length(tlist) > 0  ) {
  tblout <- rbindlist(tlist)
} else {
  write("No results found, exiting", stderr())
  quit(save = 'no', status = 0)
}
setkey(tblout, profile, accno)

# Read the annottable, if given, split into individual profiles, left join with tableout and summarise
if ( length(opt$options$annottable) > 0 ) {
  annottable <- fread(opt$options$annottable) %>% as_tibble() %>%
    # Make sure the profile column has ' & ' as separator between multiple profiles
    mutate(profile = gsub('\\s*&\\s*', ' & ', profile))
  # Check that the protein column is unique
  if ( nrow(annottable) != nrow(annottable %>% distinct(protein)) ) {
    write(sprintf("'protein' column in %s not unique, can't continue", opt$options$annottable))
    quit(save = 'no', status = 1)
  }
  # annotmap is our tool to summarise data for each entry in annottable
  annotmap   <- annottable %>% select(protein, profile) %>%
    separate_rows(profile, sep = '\\s*&\\s')
  
  # Join in the table via the annotmap and summation
  tblout <- tblout %>%
    left_join(annotmap, by = 'profile') %>% 
    group_by(accno, protein) %>% 
    summarise(profile = paste(profile, collapse = ' & '), evalue = min(evalue), score = sum(score), .groups = 'drop') %>%
    left_join(annottable %>% select(-profile), by = c('protein')) %>%
    as.data.table() %>%
    setkey(profile, accno)
}

# If we're given a scorefile as an option, read and left join
if ( length(opt$options$scorefile) > 0 ) {
  logmsg(sprintf("Joining in scores from %s, setting minimum score for absent profiles to %5.2f", opt$options$scorefile, opt$options$minscore))
  s <- read_tsv(
    opt$options$scorefile,
    col_types = cols(.default = col_character(), seq_score = col_double(), domain_score = col_double())
  ) %>%
    filter(score_type == 'GA') %>%
    select(profile, seq_score) %>%
    data.table() %>%
    setkey(profile)
  
  tblout[s, on = 'profile', seq_score := i.seq_score] 
  tblout[is.na(seq_score)]$seq_score <- opt$options$minscore
} else {
  logmsg(sprintf("Setting minimum scores to %5.2f", opt$options$minscore))
  tblout$seq_score <- opt$options$minscore
}

# Filter for scores higher than seq_score and calculate rank
logmsg("Calculating ranks")
tblout[score >= seq_score, rank := frank(desc(score)), by = accno]

logmsg(sprintf("Writing to %s", opt$options$outfile))
if ( opt$options$only_best_scoring ) {
  write_tsv(
    tblout[score >= seq_score & rank < 2] %>% select(accno, profile, rank, evalue, score) %>% arrange(accno, rank),
    opt$options$outfile
  )
} else {
  write_tsv(
    tblout[score >= seq_score] %>% select(accno, profile, rank, evalue, score) %>% arrange(accno, rank),
    opt$options$outfile
  )
}

logmsg("Done")
