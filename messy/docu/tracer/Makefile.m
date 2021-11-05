# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/main_tracer.pdf

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
	@rm -f ../pdf/main_tracer.pdf

../pdf/main_tracer.pdf:
	@pdflatex main_tracer.tex
	@pdflatex main_tracer.tex
	@pdflatex main_tracer.tex
#	@bibtex main_tracer
#	@pdflatex main_tracer.tex
#	@pdflatex main_tracer.tex
	@mv -f main_tracer.pdf ../pdf/.

../pdf/main_tracer.pdf: main_tracer.tex
