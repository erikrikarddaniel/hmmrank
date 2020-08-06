DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: hmmrank

hmmrank: hmmrank.00 hmmrank.01 hmmrank.02 hmmrank.03 hmmrank.04 hmmrank.05 hmmrank.06

hmmrank.00:
	../R/hmmrank.r --outfile=$@.out $@.d/*.tblout
	@$(CHECK)

# Empty input file
hmmrank.01:
	../R/hmmrank.r --outfile=$@.out $@.d/*.tblout
	@rc=$$?; if [ $$rc -gt 0 ]; then echo "Return code > 0, probably an error"; exit $$rc; else echo "*** $@ OK ***"; fi

hmmrank.02:
	../R/hmmrank.r --scorefile=$@.profile_scores.tsv --minscore=400 --outfile=$@.out $@.d/*.tblout
	@$(CHECK)

hmmrank.03:
	../R/hmmrank.r --outfile=$@.out $@.d/*.tblout
	@$(CHECK)

hmmrank.04:
	../R/hmmrank.r --only_best_scoring --scorefile=$@.profile_scores.tsv --minscore=400 --outfile=$@.out $@.d/*.tblout
	@$(CHECK)

hmmrank.05:
	../R/hmmrank.r --only_best_scoring --scorefile=$@.profile_scores.tsv --qfromfname --outfile=$@.out $@.d/*.tblout
	@$(CHECK)

hmmrank.06:
	../R/hmmrank.r --annottable=$@.annottable.tsv --minscore=30 --scorefile=$@.profile_scores.tsv --qfromfname --outfile=$@.out $@.d/*.tblout
	@$(CHECK)
