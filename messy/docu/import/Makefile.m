# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/main_import.pdf
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
	@rm -f ../pdf/main_import.pdf

../pdf/main_import.pdf:
	@pdflatex main_import.tex
	@pdflatex main_import.tex
	@pdflatex main_import.tex
	@bibtex main_import
	@pdflatex main_import.tex
	@pdflatex main_import.tex
	@mv -f main_import.pdf ../pdf/.

../pdf/main_import.pdf: main_import.tex main_import.bib
