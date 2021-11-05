EXTRA_DIST=Readme-MPI.txt 						\
	readme								\
	codes3								\
	MPIOM-autotools-setup.txt					\
	README.grid-setup						\
									\
	resources/doxygen_mpimet.css					\
	resources/doxygen_mpimet_footer.html				\
	resources/doxygen_mpimet_header.html				\
	resources/img/background_header.jpg				\
	resources/img/logo_mpimet.gif					\
	resources/img/logo_mpg.gif					\
									\
	tecdoc_hamocc/HAMOCC5.1_TECHNICAL_REPORT.tex			\
	tecdoc_hamocc/HAMOCC5_TECHNICAL_REPORT.tex			\
	tecdoc_hamocc/els1-4.sty					\
	tecdoc_mpiom/appendix.tex					\
	tecdoc_mpiom/bblschem.ps					\
	tecdoc_mpiom/bblschem.eps					\
	tecdoc_mpiom/bblschem.pdf					\
	tecdoc_mpiom/c1.tex						\
	tecdoc_mpiom/c2.tex						\
	tecdoc_mpiom/c3.tex						\
	tecdoc_mpiom/c4.tex						\
	tecdoc_mpiom/c5.tex						\
	tecdoc_mpiom/c6.tex						\
	tecdoc_mpiom/c7.tex						\
	tecdoc_mpiom/cgricon_new.ps					\
	tecdoc_mpiom/cgricon_new.eps					\
	tecdoc_mpiom/cgricon_new.pdf					\
	tecdoc_mpiom/grob_grid.ps					\
	tecdoc_mpiom/grob_grid.eps					\
	tecdoc_mpiom/grob_grid.pdf					\
	tecdoc_mpiom/tecdoc.tex						\
	tecdoc_mpiom/thesis.cls						\
	tecdoc_mpiom/title.tex

html-local:
	test -d resources || $(LN_S) @srcdir@/resources .
	cd .. && @DOXYGEN@
	$(MKDIR_P) html/img
	cp resources/img/* html/img 

clean-local:
	rm -r -f html
