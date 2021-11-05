# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/howto-newsubmodel.pdf

.PHONY: clean
clean:
	@rm -f *.aux
	@rm -f *.bbl
	@rm -f *.blg
	@rm -f *.log
	@rm -f *.flc
	@rm -f *.toc
	@rm -f *.out

.PHONY: distclean
distclean: clean
	@rm -f ../pdf/howto-newsubmodel.pdf

../pdf/howto-newsubmodel.pdf:
	@pdflatex howto-newsubmodel.tex
	@pdflatex howto-newsubmodel.tex
	@pdflatex howto-newsubmodel.tex
#	@bibtex howto-newsubmodel
#	@pdflatex howto-newsubmodel.tex
#	@pdflatex howto-newsubmodel.tex
	@mv -f howto-newsubmodel.pdf ../pdf/.

../pdf/howto-newsubmodel.pdf: howto-newsubmodel.tex
