BASENAME = jval_jvpp_manual
SUFFIXES = log toc lof lot dvi aux brf blg ilg ind idx vrb snm out nav
TMPFILES := $(addprefix "$(BASENAME).", $(SUFFIXES))

all: clean pdflatex 

pdflatex:
	pdflatex  $(BASENAME).tex
	bibtex    $(BASENAME).aux
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
