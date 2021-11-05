# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/main_timer.pdf

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
	@rm -f ../pdf/main_timer.pdf

../pdf/main_timer.pdf:
	@pdflatex main_timer.tex
	@pdflatex main_timer.tex
	@pdflatex main_timer.tex
	@bibtex main_timer
	@pdflatex main_timer.tex
	@pdflatex main_timer.tex
	@mv -f main_timer.pdf ../pdf/.

../pdf/main_timer.pdf: main_timer.tex main_timer.bib
../pdf/main_timer.pdf: timer_logo.pdf timer_modules.pdf
