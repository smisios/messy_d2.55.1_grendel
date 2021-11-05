# -*- makefile-gmake -*-

.NOTPARALLEL:
SHELL           = /bin/sh
#VPATH=$(SRCDIR)/mct 
VPATH=.
# SOURCE FILES

MODULE		= oasmct

PF := oas_

SRCS_F90       := m_MCTWorld.F90		\
		  m_AttrVect.F90		\
		  m_GlobalMap.F90	        \
		  m_GlobalSegMap.F90		\
		  m_GlobalSegMapComms.F90	\
		  m_Accumulator.F90		\
		  m_SparseMatrix.F90      	\
		  m_Navigator.F90		\
		  m_AttrVectComms.F90		\
		  m_AttrVectReduce.F90		\
		  m_AccumulatorComms.F90  	\
		  m_GeneralGrid.F90		\
		  m_GeneralGridComms.F90	\
		  m_SpatialIntegral.F90		\
		  m_SpatialIntegralV.F90	\
		  m_MatAttrVectMul.F90    	\
		  m_Merge.F90			\
		  m_GlobalToLocal.F90		\
		  m_ExchangeMaps.F90		\
		  m_ConvertMaps.F90		\
		  m_SparseMatrixDecomp.F90	\
		  m_SparseMatrixToMaps.F90	\
		  m_SparseMatrixComms.F90 	\
		  m_SparseMatrixPlus.F90	\
		  m_Router.F90			\
		  m_Rearranger.F90        	\
		  m_Transfer.F90                \
                  mct_mod.F90

PFSRCS_F90_A := $(subst m_,${PF}m_,${SRCS_F90})
PFSRCS_F90   := $(subst mct_,${PF}mct_,${PFSRCS_F90_A})

OBJS_ALL	= $(PFSRCS_F90:.F90=.o)

# MACHINE AND COMPILER FLAGS

include ../Makefile.conf
libdir = ../../../../../../lib

# TARGETS

all: 	$(libdir)/lib$(MODULE).a

$(libdir)/lib$(MODULE).a: $(PFSRCS_F90) $(OBJS_ALL)
	$(RM) $@
	$(AR) $@ $(OBJS_ALL)

# ADDITIONAL FLAGS SPECIFIC FOR MCT COMPILATION

MCTFLAGS = $(INCFLAG)$(MPEUPATH)

# RULES

.SUFFIXES:
.SUFFIXES: .F90 .o

$(PF)mct_mod.F90: ../../mct/mct_mod.F90
	sed "s/\($(PF)\)*mct_/$(PF)mct_/g" $< > $@
	sed -i "s/use m_/use $(PF)m_/g" $@
	sed -i "s/$(PF)m_List/m_List/g" $@
	sed -i "s/$(PF)m_string/m_string/g" $@
	sed -i "s/$(PF)m_die/m_die/g" $@
	sed -i "s/$(PF)m_MergeSort/m_MergeSort/g" $@
	sed -i "s/$(PF)m_inpak90/m_inpak90/g" $@
	sed -i "s/$(PF)m_Permuter/m_Permuter/g" $@
	sed -i "s/use m_MCTW/use $(PF)m_MCTW/g" $@
	sed -i "s/use m_Glob/use $(PF)m_Glob/g" $@
	sed -i "s/use m_AttrV/use $(PF)m_AttrV/g" $@
	sed -i "s/use m_Accum/use $(PF)m_Accum/g" $@
	sed -i "s/use m_Gener/use $(PF)m_Gener/g" $@
	sed -i "s/use m_Spa/use $(PF)m_Spa/g" $@
	sed -i "s/use m_Nav/use $(PF)m_Nav/g" $@
	sed -i "s/use m_Con/use $(PF)m_Con/g" $@
	sed -i "s/use m_Ex/use $(PF)m_Ex/g" $@
	sed -i "s/use m_Rou/use $(PF)m_Rou/g" $@
	sed -i "s/use m_Rear/use $(PF)m_Rear/g" $@

$(PF)%.F90 :: ../../mct/%.F90
	sed "s/module m_/module $(PF)m_/g" $< > $@
	sed -i "s/use m_MCTW/use $(PF)m_MCTW/g" $@
	sed -i "s/use m_Glob/use $(PF)m_Glob/g" $@
	sed -i "s/use m_AttrV/use $(PF)m_AttrV/g" $@
	sed -i "s/use m_Accum/use $(PF)m_Accum/g" $@
	sed -i "s/use m_Gener/use $(PF)m_Gener/g" $@
	sed -i "s/use m_Spa/use $(PF)m_Spa/g" $@
	sed -i "s/use m_Nav/use $(PF)m_Nav/g" $@
	sed -i "s/use m_Con/use $(PF)m_Con/g" $@
	sed -i "s/use m_Ex/use $(PF)m_Ex/g" $@
	sed -i "s/use m_Rou/use $(PF)m_Rou/g" $@
	sed -i "s/use m_Rear/use $(PF)m_Rear/g" $@

.F90.o:
	$(FC) -c $(INCPATH) -I. $(FPPDEFS) $(FCFLAGS) $(MCTFLAGS) $<

clean:
#	${RM} *.o *.mod lib$(MODULE).a
	rm -f *.o *.mod *.F90

install: all
#	$(MKINSTALLDIRS) $(libdir) $(includedir)
#	$(INSTALL) lib$(MODULE).a -m 644 $(libdir)
#	@for modfile in *.mod; do                         \
#	  echo $(INSTALL) $$modfile -m 644 $(includedir); \
#	  $(INSTALL) $$modfile -m 644 $(includedir);      \
#	done

# DEPENDENCIES

$(OBJS_ALL): $(libdir)/libmpeu.a

$(PF)m_AttrVect.o: 
$(PF)m_Accumulator.o: $(PF)m_AttrVect.o
$(PF)m_GlobalMap.o:
$(PF)m_GlobalSegMap.o:
$(PF)m_GlobalSegMapComms.o: $(PF)m_GlobalSegMap.o
$(PF)m_Navigator.o:
$(PF)m_AttrVectComms.o: $(PF)m_AttrVect.o $(PF)m_GlobalMap.o
$(PF)m_AttrVectReduce.o: $(PF)m_AttrVect.o
$(PF)m_AccumulatorComms.o: $(PF)m_AttrVect.o $(PF)m_GlobalMap.o $(PF)m_AttrVectComms.o
$(PF)m_SparseMatrix.o: $(PF)m_AttrVect.o $(PF)m_GlobalMap.o $(PF)m_AttrVectComms.o
$(PF)m_GeneralGrid.o: $(PF)m_AttrVect.o
$(PF)m_GeneralGridComms.o: $(PF)m_AttrVect.o $(PF)m_GeneralGrid.o $(PF)m_AttrVectComms.o $(PF)m_GlobalMap.o $(PF)m_GlobalSegMap.o
$(PF)m_MatAttrVectMul.o: $(PF)m_AttrVect.o $(PF)m_SparseMatrix.o $(PF)m_GlobalMap.o $(PF)m_GlobalSegMap.o $(PF)m_SparseMatrixPlus.o $(PF)m_Rearranger.o
$(PF)m_Merge.o: $(PF)m_AttrVect.o $(PF)m_GeneralGrid.o
$(PF)m_Router.o: $(PF)m_GlobalToLocal.o $(PF)m_MCTWorld.o $(PF)m_GlobalSegMap.o $(PF)m_ExchangeMaps.o
$(PF)m_Rearranger.o: $(PF)m_Router.o $(PF)m_MCTWorld.o $(PF)m_GlobalSegMap.o $(PF)m_AttrVect.o
$(PF)m_GlobalToLocal.o: $(PF)m_GlobalSegMap.o
$(PF)m_ExchangeMaps.o: $(PF)m_GlobalMap.o $(PF)m_GlobalSegMap.o $(PF)m_MCTWorld.o $(PF)m_ConvertMaps.o
$(PF)m_ConvertMaps.o: $(PF)m_GlobalMap.o $(PF)m_GlobalSegMap.o $(PF)m_MCTWorld.o
$(PF)m_SparseMatrixDecomp.o: $(PF)m_SparseMatrix.o $(PF)m_GlobalSegMap.o
$(PF)m_SparseMatrixToMaps.o: $(PF)m_SparseMatrix.o $(PF)m_GlobalSegMap.o
$(PF)m_SparseMatrixComms.o:	$(PF)m_SparseMatrix.o $(PF)m_SparseMatrixDecomp.o $(PF)m_GlobalSegMap.o $(PF)m_AttrVectComms.o
accumulate.o: $(PF)m_AttrVect.o $(PF)m_Accumulator.o
$(PF)m_SpatialIntegral.o: $(PF)m_SpatialIntegralV.o $(PF)m_GeneralGrid.o $(PF)m_AttrVect.o $(PF)m_AttrVectReduce.o
$(PF)m_SpatialIntegralV.o: $(PF)m_AttrVect.o $(PF)m_AttrVectReduce.o
$(PF)m_Transfer.o: $(PF)m_AttrVect.o $(PF)m_Router.o $(PF)m_MCTWorld.o
$(PF)m_SparseMatrixPlus.o: $(PF)m_GlobalSegMap.o $(PF)m_Rearranger.o $(PF)m_SparseMatrix.o $(PF)m_SparseMatrixComms.o $(PF)m_SparseMatrixToMaps.o $(PF)m_GlobalToLocal.o
$(PF)mct_mod.o:  $(PF)m_Accumulator.o  $(PF)m_AttrVect.o  $(PF)m_AttrVectComms.o $(PF)m_GeneralGrid.o $(PF)m_GeneralGridComms.o $(PF)m_GlobalSegMap.o $(PF)m_GlobalSegMapComms.o $(PF)m_GlobalToLocal.o $(PF)m_MatAttrVectMul.o $(PF)m_MCTWorld.o $(PF)m_Rearranger.o $(PF)m_Router.o $(PF)m_SparseMatrix.o $(PF)m_SparseMatrixComms.o $(PF)m_SparseMatrixPlus.o $(PF)m_SparseMatrixToMaps.o $(PF)m_Transfer.o
