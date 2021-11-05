# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/main_grid.pdf

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
	@rm -f ../pdf/main_grid.pdf

../pdf/main_grid.pdf:
	@pdflatex main_grid.tex
	@pdflatex main_grid.tex
	@pdflatex main_grid.tex
	@bibtex main_grid
	@pdflatex main_grid.tex
	@pdflatex main_grid.tex
	@mv -f main_grid.pdf ../pdf/.

../pdf/main_grid.pdf: main_grid.tex main_grid.bib
../pdf/main_grid.pdf: grid_modules.pdf
