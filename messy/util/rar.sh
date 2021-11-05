#!/bin/sh -e

# mz_sg_20170322(@tonnerre.mpic.de)
# info: this utility script creates a rar archive of the distribution
# based on: zip.sh
# output file format: YYMMDD.<distro_name>.rar

# check if RAR is available
if ! type "rar" > /dev/null; then
   echo "Cannot detect RAR installed on the system. Will not continue."
   exit 1
fi

# if the current directory is only a link, the following cd command
# will jump into the real directory:
cd `/bin/pwd`

# DO NOT RAR VECTORIZED SMIL FILE
VEC=0
if test -L  ./messy/echam5/smil/messy_mecca1_e5.f90 ; then
   echo '... resetting messy/echam5/smil/messy_mecca1_e5.f90 for raring'
   mv -f ./messy/echam5/smil/messy_mecca1_e5.f90 \
         ./messy/echam5/smil/messy_mecca1_e5.f90.old
   mv -f ./messy/echam5/smil/messy_mecca1_e5.f90-ori \
         ./messy/echam5/smil/messy_mecca1_e5.f90
   VEC=1
fi


DIR=`pwd | sed 's|/.*/||'`
FIL=`date +%y%m%d`.$DIR'_src'
if test -r $FIL.rar ; then
  # mz_rs_20100417+
  #echo "rar file exists! Please (re)move it first:"
  #echo "rm $FIL.rar" 
  #exit 1
  echo "Renaming old $FIL.rar file to $FIL.rar~"
  mv -f $FIL.rar $FIL.rar~
  # mz_rs_20100417-
fi
cd ..
if test -r $FIL.rar ; then
  # mz_rs_20100417+
  #echo "rar file exists! Please (re)move it first:"
  #echo "rm ../$FIL.rar" 
  #exit 1
  echo "Renaming old $FIL.rar file to $FIL.rar~"
  mv -f $FIL.rar $FIL.rar~
  # mz_rs_20100417-
fi

# rar options:
# -tl         set archive time to latest file
# -r          recurse subdirectories
# -ol[a]/-oh  process symbolic/hard links as the link
# -rr         add recovery record
# -x'...'     exclude files

ARC='rar a -tl -r -ol -rr'

$ARC \
 -x'*~' -x'*.mod' -x'*.exe' -x'*.o' -x'*.a' -x'*.nc' -x'*.md5' \
 -x'*.log' -x'*.old' -x'*/ferret.jnl' -x'*.zip' -x'*.rar' -x'*.tar' \
 -x'*/i.*.L' -x'*/i.*.O' -x'*/F*.f' -x'*/*.lst' -x'*/F*.f90' -x'*.dat' \
 -x'*/fort.*' -x'*/tmp_*' -x'*/*.i90' -x'*/*.i' -x'*/*.optrpt'\
 -x'*/*.d' \
 -x'*.pyc' -x'*.aux' -x'*.bbl' -x'*.toc' -x'*.blg' -x'*/.git/*' \
 -x'*/config/config.h' -x'*/config.cache' -x'*/config.status' \
 -x$DIR/'Makefile' -x$DIR/'intera*/Makefile' -x$DIR/'echam5*/Makefile' \
 -x$DIR/'mpiom*/Makefile' -x$DIR/'cosmo*/Makefile' \
 -x$DIR/'cesm1*/Makefile' \
 -x$DIR/'libsrc/pio/Makefile.conf' -x$DIR/'libsrc/mct/Makefile.conf' \
 -x'*/workdir*/?*' \
 -x$DIR/'echam5*/__*' \
 -x$DIR/'cesm1*/__*' \
 -x$DIR/'mpiom*/__*' \
 -x$DIR/'cosmo*/__*' \
 -x$DIR/'messy/util/*.pl' \
 -x$DIR/'messy/mbm/*/__*' \
 -x$DIR/'messy/tools/kpp/bin/?*' \
 -x$DIR/'messy/tools/kpp1/bin/?*' \
 -x$DIR/'messy/util/locate_f90.sh' \
 -x$DIR/'messy/echam5/smil/*mecca1_vec*' \
 -x$DIR/'messy/smcl/*mecca1_vec*' \
 -x$DIR/'messy/mbm/caaba/output/?*' \
 -x$DIR/'messy/mbm/caaba/skeleton/output/?*' \
 -x$DIR/'messy/mbm/mecca1/vector/smcl/*kpp*' \
 -x$DIR/'messy/mbm/mecca1/vector/smil/*.inc' \
 -x$DIR/'messy/tools/kpp1/bin/kpp' \
 -x$DIR/'messy/tools/kpp/bin/kpp' \
 -x'*/depend.mk' \
 -x'*/messy/bin/?*' \
 -x'*/tmp/*' \
 -x$DIR/'messy/smcl/messy_main_compilerinfo_mem.f90' \
 -x$DIR/'messy/tools/jvpp/dat_lit/old/workdir_176' \
 -x$DIR/'messy/tools/jvpp/dat_m17/old/workdir_176' \
$FIL.rar $DIR

#  mz_rs_20100727+
# add input files for CAABA (even though they are *.nc files):
$ARC $FIL.rar $DIR/messy/mbm/caaba/traject/*.nc
$ARC $FIL.rar $DIR/messy/mbm/caaba/multirun/input/*.nc
$ARC $FIL.rar $DIR/messy/mbm/caaba/input/*.nc
#  mz_rs_20100727-

#  op_pj_20161108+
if test -d $DIR/icon* ; then
  $ARC $FIL.rar $DIR/icon*/data/*
fi
#  op_pj_20161108-

mv -f $FIL.rar $DIR/.

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
