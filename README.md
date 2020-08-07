# hmmrank
Calculates a ranking of hmm profiles per sequence

Given a set of hmmsearch "tblout" result files, this R script outputs a table where each sequence hits
against a set of hmm profiles are ranked based on score. It requires the `optparse` package as well as
`readr`, `tidyr` and `dplyr` from Tidyverse.

## Installation and documentation

After installing the package with Conda, you can run the script with`--help` to 
get help...

```
$ conda create -n hmmrank -c matricaria hmmrank
$ hmmrank.r --help
Usage: ../R/hmmrank.r [options] file0.tblout .. filen.tblout

  Ranks tabular output from hmmsearch (HMMER package). Optionally filters based on score and adds annotation information.

  This is a brief overview of options, check out https://github.com/erikrikarddaniel/hmmrank/blob/master/README.md for more instructions.


Options:
        --annottable=ANNOTTABLE
                Name of table with annotation assignments, at a minium must contain 'protein' and 'profile' columns. The 'profile' column might contain multiple profiles, separated by '&'.

        --minscore=MINSCORE
                Minimum score, default 0

        --only_best_scoring
                Output only the best scoring, i.e. scores < 2

        --outfile=OUTFILE
                Name of output file

        --qfromfname
                Extract query name from file name -- basename until first period -- instead of taking it from the query column in each data file, default false.

        --scorefile=SCOREFILE
                Name of tsv file containing accept scores (GA in score_type): profile, score_type and seq_score

        -v, --verbose
                Print progress messages

        -V, --version
                Print script version number

        -h, --help
                Show this help message and exit

```

### Running `hmmrank` with an annotation table


#### Entries with multiple profiles

The `profile` column may contain entries like `TIGR01466 & TIGR01467` which will be taken as an
entry that requires both profiles to match a sequence. Scores for the two matches will be summed and
the combined entry will hence receive a better score than any of the two by themselves and the
combination will be ranked higher. 

Note that combinations that are observed in the data but not mentioned in the annotation table will
be disregarded. E.g., if `PF00590` and `TIGR01469` both appear in the annotation table, proteins
that match both profiles will be assigned to whichever of the two they match best.

Also note that score filtering, if asked for, takes place before evaluating combinations. Long
proteins that match multiple profiles may be penalized by this if the target sequence is a fragment.
You may want to lower high minimum scores to deal with this problem; take a particular look at
TIGRFAMs since they typically are long proteins with high default GA scores.
