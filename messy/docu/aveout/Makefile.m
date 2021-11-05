# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/aveout.pdf

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
	@rm -f ../pdf/aveout.pdf

../pdf/aveout.pdf:
	@pdflatex aveout.tex
	@pdflatex aveout.tex
	@pdflatex aveout.tex
	@mv -f aveout.pdf ../pdf/.

../pdf/aveout.pdf: aveout.tex

