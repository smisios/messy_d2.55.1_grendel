#!/bin/sh -e

cd ..
DIR=`pwd | sed 's|/.*/||'`
FIL=$DIR'_src'
if test -r $FIL.zip ; then
  echo "Renaming old $FIL.zip file to $FIL.zip~"
  mv -f $FIL.zip $FIL.zip~
fi
cd ..
if test -r $FIL.zip ; then
  echo "Renaming old $FIL.zip file to $FIL.zip~"
  mv -f $FIL.zip $FIL.zip~
  # mz_rs_20100417-
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
 -x '*/fort.*' -x '*/tmp_*' -x '*/*.i90' -x '*/*.i' -x '*/*.optrpt' \
 -x '*/*.d' \
 -x '*core.?????' \
 -x '*.pyc' -x '*.aux' -x '*.bbl' -x '*.toc' -x '*.blg' -x '*/.git/*' \
 -x '*/config/config.h' -x '*/config.cache' -x '*/config.status' \
 -x $DIR/'Makefile' -x '*/workdir*/?*' \
 -x $DIR/src/'cclm*/__*' \
 -x $DIR/'data/gcm*' \
 -x $DIR/'data/ext*' \
 -x $DIR/'chain/gcm_to_cclm/sp???/jobs/*.job' \
 -x $DIR/'chain/gcm_to_cclm/sp???/jobs/functions.sh' \
 -x $DIR/'chain/cclm_to_cclm/sp???/jobs/*.job*' \
 -x $DIR/'chain/cclm_to_cclm/sp???/jobs/functions.sh' \
 -x $DIR/'chain/gcm_to_cclm/sp???/utils/samoa*' \
 -x $DIR/'chain/cclm_to_cclm/sp???/utils/samoa*' \
 -x $DIR/'chain/arch/sp???/*' \
 -x $DIR/'chain/work/sp???/*' \
 -x $DIR/'chain/scratch/sp???/*' \
 -x $DIR/'step_by_step/gcm_to_cclm/sp_*.o*' \
 -x $DIR/'step_by_step/gcm_to_cclm/log_cclm/*' \
 -x $DIR/'step_by_step/gcm_to_cclm/log_int2lm/*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/bin*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/import*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/nml*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/save' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/*.nml' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/M*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/YU*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/INPUT*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/xmessy*' \
 -x $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/restarts/*' \
 -x $DIR/'step_by_step/cclm_to_cclm/sp_*.o*' \
 -x $DIR/'step_by_step/cclm_to_cclm/log_cclm/*' \
 -x $DIR/'step_by_step/cclm_to_cclm/log_int2lm/*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/bin*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/import*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/nml*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/save' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/*.nml' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/M*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/YU*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/INPUT*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/xmessy*' \
 -x $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/restarts/*' \
 -x $DIR/src/'messy/util/*.pl' \
 -x $DIR/src/'messy/mbm/*/__*' \
 -x $DIR/src/'messy/tools/kpp/bin/?*' \
 -x $DIR/src/'messy/tools/kpp1/bin/?*' \
 -x $DIR/src/'messy/util/locate_f90.sh' \
 -x $DIR/src/'messy/echam5/smil/*mecca1_vec*' \
 -x $DIR/src/'messy/smcl/*mecca1_vec*' \
 -x $DIR/src/'messy/mbm/caaba/output/?*' \
 -x $DIR/src/'messy/mbm/caaba/skeleton/output/?*' \
 -x $DIR/src/'messy/mbm/mecca1/vector/smcl/*kpp*' \
 -x $DIR/src/'messy/mbm/mecca1/vector/smil/*.inc' \
 -x $DIR/src/'messy/tools/kpp1/bin/kpp' \
 -x $DIR/src/'messy/tools/kpp/bin/kpp' \
 -x '*/depend.mk' \
 -x '*/messy/bin/?*' \
 -x '*/tmp/*' \
 -x $DIR/src/'messy/smcl/messy_main_compilerinfo_mem.f90' \
 -x $DIR/src/'messy/tools/jvpp/dat_lit/old/workdir_176' \
 -x $DIR/src/'messy/tools/jvpp/dat_m17/old/workdir_176'


# add empty directories and .gitignore files:
zip -o $FIL.zip $DIR/'step_by_step/gcm_to_cclm/log_int2lm/.gitignore'
zip -o $FIL.zip $DIR/'step_by_step/gcm_to_cclm/log_cclm/.gitignore'
zip -o $FIL.zip $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/restarts/.gitignore'
zip -o $FIL.zip $DIR/'step_by_step/cclm_to_cclm/log_int2lm/.gitignore'
zip -o $FIL.zip $DIR/'step_by_step/cclm_to_cclm/log_cclm/.gitignore'
zip -o $FIL.zip $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/restarts/.gitignore'
zip -o $FIL.zip $DIR/'step_by_step/gcm_to_cclm/log_int2lm/'
zip -o $FIL.zip $DIR/'step_by_step/gcm_to_cclm/log_cclm/'
zip -o $FIL.zip $DIR/'step_by_step/gcm_to_cclm/data/cclm_output/restarts/'
zip -o $FIL.zip $DIR/'step_by_step/cclm_to_cclm/log_int2lm/'
zip -o $FIL.zip $DIR/'step_by_step/cclm_to_cclm/log_cclm/'
zip -o $FIL.zip $DIR/'step_by_step/cclm_to_cclm/data/cclm_output/restarts/'

mv -f $FIL.zip $DIR/.

exit 0
