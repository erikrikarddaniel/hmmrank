# hmmrank
Calculates a ranking of hmm profiles per sequence

Given a set of hmmsearch "tblout" result files, this R script outputs a table where each sequence hits
against a set of hmm profiles are ranked based on score. It requires the `optparse` package as well as
`readr`, `tidyr` and `dplyr` from Tidyverse.

## Installation and documentation

After cloning the repository and installing required libraries, you can run the script with`--help` to 
get help...

```
$ Rscript -e "install.packages(c('optparse', 'readr', 'tidyr', 'dplyr'))"
$ git clone https://github.com/erikrikarddaniel/hmmrank.git
$ hmmrank/src/R/hmmrank.r --help
Usage: hmmrank/src/R/hmmrank.r [options] file0.tblout .. filen.tblout


Options:
        --minscore=MINSCORE
                Minimum score, default 0

        --outfile=OUTFILE
                Name of output file

        --scorefile=SCOREFILE
                Name of tsv file containing accept scores (GA in score_type): profile, score_type and seq_score

        -v, --verbose
                Print progress messages

        -h, --help
                Show this help message and exit

```
