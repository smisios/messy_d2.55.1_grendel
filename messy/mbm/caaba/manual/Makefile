BASENAME = caaba_mecca_manual
SUFFIXES = log toc lof lot dvi aux brf blg ilg ind idx vrb snm out nav
TMPFILES := $(addprefix "$(BASENAME).", $(SUFFIXES))

all: clean pdflatex 

pdflatex:
	pdflatex  $(BASENAME).tex
	bibtex    $(BASENAME).aux
	pdflatex  $(BASENAME).tex
#	-s uses index style file
#       -c compresses double blanks and cr's in idx file
	makeindex $(BASENAME).idx -s myindex.ist -c
	pdflatex  $(BASENAME).tex
	pdflatex  $(BASENAME).tex

.PHONY: clean
clean:
	@for TMPFILE in $(TMPFILES) ;\
	  do \
	    if test -r $$TMPFILE ; then \
	      echo "removing temporary file $$TMPFILE" ; \
	      rm -f $$TMPFILE ; \
	    fi ; \
	  done 
