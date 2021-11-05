# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

.PHONY: docu
docu: ../pdf/MMDUserManual.pdf

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
	@rm -f ../pdf/MMDUserManual.pdf

../pdf/MMDUserManual.pdf:
	@pdflatex MMDUserManual.tex
	@pdflatex MMDUserManual.tex
	@pdflatex MMDUserManual.tex
	@bibtex MMDUserManual
	@pdflatex MMDUserManual.tex
	@pdflatex MMDUserManual.tex
	@mv -f MMDUserManual.pdf ../pdf/.

../pdf/MMDUserManual.pdf: MMDUserManual.tex MMDUM.bib
../pdf/MMDUserManual.pdf: MECOn_coupling_stecker.pdf  MMDUM_ext_i2c_c_grids.pdf  MMDUM_flowchart_iniphase.pdf  MMDUM_halo_a.pdf  MMDUM_halo_c.pdf   MMDUM_ptr.pdf MMDUM_decisiontree.pdf  MMDUM_fieldflow.pdf  MMDUM_flowchart_timeloop.pdf  MMDUM_halo_b.pdf  MMDUM_ixxxcos.pdf
