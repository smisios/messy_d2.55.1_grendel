# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

INSTALLDIR = ../../../bin/.
CHEMPROPDIR = ../../mbm/caaba/mecca/tracer/chemprop

FPC_EXE := $(shell which fpc 2>/dev/null || echo "<none>")
FPC_VER := $(shell [ "$(FPC_EXE)" != "<none>" ] && $(FPC_EXE) -iV || echo "<none>")
FPC_VER_MAJ := $(shell echo "$(FPC_VER)" | awk -F. '{print $$1}' )
FPC_VER_MIN := $(shell echo "$(FPC_VER)" | awk -F. '{$$1=""; print $$0}' )

# ----------------------------------------------
### name of the executable that will be produced
# OBJECT_MODE=32 is a fix for fpc @AIX systems
FPC   = env OBJECT_MODE=32 $(FPC_EXE) -l -viwnh -B
PROG1 = imtag.exe
PROG2 = embudget.exe
CFG   = $(shell test -d $(CHEMPROPDIR) && echo TRACDEF_CHEMPROP || echo TRACDEF_ADD_H2O )
# the option below enables additional diagnostic for unaccounted mass prod./loss
#CFG  += USE_PT_UPL

# --------------------------------------------------------------------
DEFOPT = -d
CDEFS = $(addprefix $(DEFOPT), $(CFG))

all: $(PROG1) $(PROG2)
	@[ "$(FPC_EXE)" != "<none>" ] && echo "CFG: $(CFG)" || /bin/true

$(PROG1): imtag.pas imcom.inc imdot.inc checkfpc
	@if [ "$(FPC_EXE)" != "<none>" ]; then \
	  $(FPC) $(CDEFS) imtag -o$(PROG1) ;\
	  mv -f $(PROG1) $(INSTALLDIR)/. ;\
	fi

$(PROG2): embudget.pas imcom.inc checkfpc
	@if [ "$(FPC_EXE)" != "<none>" ]; then \
	  $(FPC) $(CDEFS) embudget -o$(PROG2) ;\
	  mv -f $(PROG2) $(INSTALLDIR)/. ;\
	fi

ctags:
	@rm -f TAGS tags
	ctags --language-force=pascal --fields=+iaS --extra=+q -e *.pas *.inc

debug:
	$(FPC) $(CDEFS) -g -gh -gl imtag -o$(PROG1)
	mv -f $(PROG1) $(INSTALLDIR)/.
	$(FPC) $(CDEFS) -g -gh -gl embudget -o$(PROG2)
	mv -f $(PROG2) $(INSTALLDIR)/.

checkfpc:
	@if [ "$(FPC_EXE)" = "<none>" ]; then \
	  echo '###########################################################';\
	  echo '### Warning: fpc (free pascal compiler) not available.';\
	  echo '###########################################################';\
	  exit 0 ;\
	else \
	echo "fpc version@path: $(FPC_VER)@$(FPC_EXE)"; \
	  if [ $(FPC_VER_MAJ) -lt 3 ]; then \
	    echo '###########################################################';\
	    echo '### Error: fpc (free pascal compiler) version is too old ';\
	    echo '### Minimum required version is 3, detected is: $(FPC_VER_MAJ) ($(FPC_VER))';\
	    echo '###########################################################';\
	    exit 3 ;\
	  fi \
	fi

clean:
	rm -f *~
	rm -f *.o

distclean: clean
	rm -f $(INSTALLDIR)/$(PROG1)
	rm -f $(INSTALLDIR)/$(PROG2)

install: all

# ----------------------------------------------
