# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/main_channel.pdf

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
	@rm -f ../pdf/main_channel.pdf

../pdf/main_channel.pdf:
	@pdflatex main_channel.tex
	@pdflatex main_channel.tex
	@pdflatex main_channel.tex
	@bibtex main_channel
	@pdflatex main_channel.tex
	@pdflatex main_channel.tex
	@mv -f main_channel.pdf ../pdf/.

../pdf/main_channel.pdf: main_channel.tex main_channel.bib
../pdf/main_channel.pdf: channel_logo.pdf channel_modules.pdf
