DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: hmmrank.00

hmmrank.00:
	../R/hmmrank.r --verbose --outfile=$@.out $@.d/*.tblout
	@$(CHECK)
