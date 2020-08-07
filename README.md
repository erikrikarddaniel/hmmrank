# hmmrank

Calculates a ranking of hmm profiles per sequence.

Given a set of hmmsearch "tblout" result files, this R script outputs a table where each sequence hits
against a set of hmm profiles are ranked based on score.

## Installation and documentation

After installing the package with Conda, you can run the script with`--help` to 
get help...

```
$ conda install -c matricaria hmmrank # (Better in an environment, see Conda documentation)
$ hmmrank.r --help
Usage: ../R/hmmrank.r [options] file0.tblout .. filen.tblout

  Ranks tabular output from hmmsearch (HMMER package). Optionally filters based on score and adds annotation information.

  This is a brief overview of options, check out https://github.com/erikrikarddaniel/hmmrank/blob/master/README.md for more instructions.


Options:
        --annottable=ANNOTTABLE
                Name of table with annotation assignments, at a minium must contain 'protein' and 'profile' columns. The 'profile' column might contain multiple profiles, separated by '&'.

        --maxscore=MAXSCORE
                Maximum sequence score; overrides the sequence score in scorefile if higher than this value, default Inf

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

## Creating a scorefile

A file for the `--scorefile` parameter is a tab separated file looking like this:

```
profile	score_type	seq_score	domain_score
PF00590	GA	27.80	27.80
PF06180	GA	26.20	26.20
TIGR01466	GA	100.00	100.00
TIGR01467	GA	100.00	100.00
TIGR01469	GA	100.00	100.00
```

It can be created from a set of Pfam and TIGRFAM hmm files like this:

```
echo -e "profile\tscore_type\tseq_score\tdomain_score" > scores.tsv; grep '^GA' *.hmm | sed 's/.hmm:/\t/' | sed 's/ \+/\t/' | sed 's/;//' >> scores.tsv
```

(If there are other hmm files in the directory, that do not contain GA scores, it doesn't hurt, but
they will not be included in the output.)

The above example and command implies that you use the `--qfromfname` option to the script and that
your output files are named after the profile, minus the `.hmm` and plus `.tblout`, e.g.
`PF00590.tblout`

The score file doesn't have to include all profiles, any that are not there will be set to 0 or
whatever you specify with the `--minscore` parameter to the script.

If you have a lot of fragmentary sequences, scores in hmm files can be a little too strict (e.g.
TIGRFAMs). You can use `--maxscore` to limit the minimum score globally. In the output, the correct
score will appear so you can filter afterwards.

## Running `hmmrank` with an annotation table

You can provide a table with annotation data using the `--annottable` option. The table must contain
a `protein` and a `profile` column (names are case sensitive), but can contain any number of other
columns, which will be added to the output table.

### Entries with multiple profiles

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
TIGRFAMs since they typically are long proteins with high default GA scores. See also description of
`--maxscore` above.
