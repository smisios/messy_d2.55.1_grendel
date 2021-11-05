# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/ncregrid.pdf

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
	@rm -f ../pdf/ncregrid.pdf

../pdf/ncregrid.pdf:
	@pdflatex ncregrid.tex
	@pdflatex ncregrid.tex
	@pdflatex ncregrid.tex
#	@bibtex ncregrid
#	@pdflatex ncregrid.tex
#	@pdflatex ncregrid.tex
	@mv -f ncregrid.pdf ../pdf/.

../pdf/ncregrid.pdf: ncregrid.tex
