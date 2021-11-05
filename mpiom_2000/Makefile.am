SUBDIRS=support bin doc

if SETUPPROGS
SUBDIRS += setup
endif

ACLOCAL_AMFLAGS=-I m4

EXTRA_DIST=Doxyfile-1.5.6					\
	Doxyfile-1.5.8						\
	DoxygenLayout.xml					\
	.gitignore						\
	INSTALLATION						\
	LICENSE							\
	PROBLEMS						\
								\
	config/mh-aix 						\
	config/mh-alpha						\
	config/mh-cray-x1					\
	config/mh-cray-xd1					\
	config/mh-cray-xt3					\
	config/mh-darwin					\
	config/mh-es						\
	config/mh-hpux						\
	config/mh-irix6						\
	config/mh-linux						\
	config/mh-linux-ia64					\
	config/mh-linux-x64 					\
	config/mh-nec						\
	config/mh-solaris					\
								\
	contrib/aix/README					\
								\
	plots/2dplot/gr2ps_bash					\
	plots/2dplot/mpiom.partab				\
	plots/timeser/documentation.pdf				\
	plots/timeser/documentation.tex				\
	plots/timeser/tex2pdf					\
	plots/timeser/zeitser.partab				\
	plots/timeser/tsplot					\
								\
	run/prepare_sge_run_mpiom_test				\
	run/prepare_sge_run_mpiom_core				\
								\
	util/createMakefiles.pl.in				\
	util/doxyfy						\
	util/expand_template					\
	util/fgenrc						\
	util/foreach.m4						\
	util/sunf95preproc-wrapper				\
	util/svnfileversion					\
	util/sxpreproc-wrapper					\
	util/xlfpreproc-wrapper

doc: html

index:
	\rm -f html/[a-z]*
	util/f2html.pl -f util/fgenrc -d html @top_srcdir@/src_hamocc \
	@top_srcdir@/src

clean-local:
	\rm -f html/[a-z]*
