# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/scalc.pdf

.PHONY: clean
clean:
	@rm -f *.aux
	@rm -f *.bbl
	@rm -f *.blg
	@rm -f *.log
	@rm -f *.flc
	@rm -f *.toc

.PHONY: distclean
distclean: clean
	@rm -f ../pdf/scalc.pdf

../pdf/scalc.pdf:
	@pdflatex scalc.tex
	@pdflatex scalc.tex
	@pdflatex scalc.tex
	@bibtex scalc
	@pdflatex scalc.tex
	@pdflatex scalc.tex
	@mv -f scalc.pdf ../pdf/.

../pdf/scalc.pdf: scalc.tex scalc.bib scalc_front.pdf
