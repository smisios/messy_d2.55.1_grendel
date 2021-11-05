#!/bin/sh -e

# if the current directory is only a link, the following cd command
# will jump into the real directory:
cd `/bin/pwd`

# DO NOT ZIP VECTORIZED SMIL FILE
VEC=0
if test -L  ./messy/echam5/smil/messy_mecca1_e5.f90 ; then
   echo '... resetting messy/echam5/smil/messy_mecca1_e5.f90 for zipping'
   mv -f ./messy/echam5/smil/messy_mecca1_e5.f90 \
	 ./messy/echam5/smil/messy_mecca1_e5.f90.old
   mv -f ./messy/echam5/smil/messy_mecca1_e5.f90-ori \
	 ./messy/echam5/smil/messy_mecca1_e5.f90
   VEC=1
fi

ODIR=`pwd | sed 's|/.*/||'`

if test -d .git ; then
    FIL='messy_'
    FIL=${FIL}`git describe --abbrev=0`
    FIL=${FIL}_`git branch | grep \* | cut -d ' ' -f2`
    FIL=${FIL}_`git rev-parse --verify --quiet HEAD | awk '{print substr($0,1,8)}'`
    set +e
    MOD="`git status --porcelain | grep -E '(f90|inc)$$'`"
    set -e
    if test "$MOD" != "" ; then
      FIL=${FIL}'_mod'
    fi
    DIR=$FIL
    FIL=$FIL'_src'
else
    DIR=$ODIR
    FIL=$DIR'_src'
fi
#echo $ODIR
#echo $DIR
#echo $FIL

if test -r $FIL.zip ; then
  # mz_rs_20100417+
  #echo "zip file exists! Please (re)move it first:"
  #echo "rm $FIL.zip"
  #exit 1
  echo "Renaming old $FIL.zip file to $FIL.zip~"
  mv -f $FIL.zip $FIL.zip~
  # mz_rs_20100417-
fi
cd ..
if test -r $FIL.zip ; then
  # mz_rs_20100417+
  #echo "zip file exists! Please (re)move it first:"
  #echo "rm ../$FIL.zip"
  #exit 1
  echo "Renaming old $FIL.zip file to $FIL.zip~"
  mv -f $FIL.zip $FIL.zip~
  # mz_rs_20100417-
fi

if test "$DIR" != "$ODIR" ; then
   if test -e $DIR ; then
      echo "directory exists! Please (re)move it first: ../"$DIR
      exit 1
   fi
   mv $ODIR $DIR
fi

# zip options:
# -o make zipfile as old as latest entry
# -r recurse into directories
# -y store symbolic links as the link instead of the referenced file
# -x '...' exclude files

zip -ory $FIL.zip $DIR \
 -x '*~' -x '*.mod' -x '*.exe' -x '*.o' -x '*.a' -x '*.nc' -x '*.md5' \
 -x '*.log' -x '*.old' -x '*/ferret.jnl' -x '*.zip' -x '*.tar' \
 -x '*/i.*.L' -x '*/i.*.O' -x '*/F*.f' -x '*/*.lst' -x '*/F*.f90' -x '*.dat' \
 -x '*/fort.[0-9]*' -x '*/tmp_*' -x '*/*.i90' -x '*/*.i' -x '*/*.optrpt' \
 -x '*/*.d' \
 -x '*.pyc' -x '*.aux' -x '*.bbl' -x '*.toc' -x '*.blg' -x '*/.git/*' \
 -x '*/.git*' \
 -x '*/config/config.h' -x '*/config.cache' -x '*/config.status' \
 -x $DIR/'Makefile' -x $DIR/'intera*/Makefile' -x $DIR/'echam5*/Makefile' \
 -x $DIR/'mpiom*/Makefile' -x $DIR/'cosmo*/Makefile' \
 -x $DIR/'cesm1*/Makefile' \
 -x $DIR/'libsrc/pio/Makefile.conf' -x $DIR/'libsrc/mct/Makefile.conf' \
 -x '*/workdir*/?*' \
 -x $DIR/'echam5*/__*' \
 -x $DIR/'cesm1*/__*' \
 -x $DIR/'mpiom*/__*' \
 -x $DIR/'cosmo*/__*' \
 -x $DIR/'messy/util/*.pl' \
 -x $DIR/'messy/mbm/*/__*' \
 -x $DIR/'messy/tools/kpp/bin/?*' \
 -x $DIR/'messy/tools/kpp1/bin/?*' \
 -x $DIR/'messy/util/locate_f90.sh' \
 -x $DIR/'messy/echam5/smil/*mecca1_vec*' \
 -x $DIR/'messy/smcl/*mecca1_vec*' \
 -x $DIR/'messy/mbm/caaba/output/?*' \
 -x $DIR/'messy/mbm/caaba/skeleton/output/?*' \
 -x $DIR/'messy/mbm/mecca1/vector/smcl/*kpp*' \
 -x $DIR/'messy/mbm/mecca1/vector/smil/*.inc' \
 -x $DIR/'messy/tools/kpp1/bin/kpp' \
 -x $DIR/'messy/tools/kpp/bin/kpp' \
 -x '*/depend.mk' \
 -x '*/messy/bin/?*' \
 -x '*/tmp/*' \
 -x $DIR/'messy/smcl/messy_main_compilerinfo_mem.f90' \
 -x $DIR/'messy/tools/jvpp/dat_lit/old/workdir_176' \
 -x $DIR/'messy/tools/jvpp/dat_m17/old/workdir_176' \
 -x $DIR/'libsrc/oasis3-mct/lib/mct/build' \
 -x $DIR/'libsrc/oasis3-mct/lib/psmile/build'

# special files added in otherwise ignored directories
if test -r libsrc/oasis3-mct ; then
zip -oy $FILE.zip $DIR/'libsrc/oasis3-mct/lib/psmile/build/Makefile.m'
fi

#  mz_rs_20180320+
# add several files for CAABA that have been skipped so far:
zip -oy $FIL.zip $DIR/messy/mbm/caaba/skeleton/samplepoints/*.nc
zip -oy $FIL.zip $DIR/messy/mbm/caaba/mecca/graphtool/caaba_mecca_rr.nc
zip -oy $FIL.zip $DIR/messy/mbm/caaba/pycaaba/__init__.py
#  mz_rs_20180320-

#  op_pj_20161108+
#  zip -o $FIL.zip $DIR/icon*/data/*
idirs="`find $DIR -maxdepth 1 -type d -name 'icon*'`"
if test "$idirs" != "" ; then
   zip -o $FIL.zip $DIR/icon*/data/*
fi
#  op_pj_20161108-

mv -f $FIL.zip $DIR/.

if test "$DIR" != "$ODIR" ; then
   mv $DIR $ODIR
fi

# RESET LOCAL VEC-STATUS
if test "$VEC" = "1" ; then
   echo '... resetting messy/echam5/smil/messy_mecca1_e5.f90'
   cd -
   mv -f ./messy/echam5/smil/messy_mecca1_e5.f90 \
	 ./messy/echam5/smil/messy_mecca1_e5.f90-ori
   mv -f ./messy/echam5/smil/messy_mecca1_e5.f90.old \
	 ./messy/echam5/smil/messy_mecca1_e5.f90
fi

exit 0
