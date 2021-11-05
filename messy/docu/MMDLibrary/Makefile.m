# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/MMDLibrary.pdf

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
	@rm -f ../pdf/MMDLibrary.pdf

../pdf/MMDLibrary.pdf:
	@pdflatex MMDLibrary.tex
	@pdflatex MMDLibrary.tex
	@pdflatex MMDLibrary.tex
	@bibtex MMDLibrary
	@pdflatex MMDLibrary.tex
	@pdflatex MMDLibrary.tex
	@mv -f MMDLibrary.pdf ../pdf/.

../pdf/MMDLibrary.pdf: MMDLibrary.tex MMDlib.bib
../pdf/MMDLibrary.pdf: MMDlib_files.pdf MMDlib_idxdef_server.pdf  MMDlib_MPIcomm.pdf  MMDlib_SC_hierarchy.pdf MMDlib_idxdef_client.pdf  MMDlib_model_layout.pdf   MMDlib_NrEle.pdf  MMDlib_work_flow.pdf
