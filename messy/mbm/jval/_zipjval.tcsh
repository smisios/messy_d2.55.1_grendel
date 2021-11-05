#! /bin/tcsh -f
# Time-stamp: <2014-01-15 16:01:41 sander>
# _zipjval.tcsh: create a zip file of jval code

# if the current directory is only a link, the following cd command
# will jump into the real directory:
cd `/bin/pwd`

if ( "$1" == "" ) then
  echo "This script should be called via the Makefile, e.g.:"
  echo "  gmake zip    --> archive important files"
  echo "  gmake zipall --> archive all files"
  exit
endif

set dirname = `basename $PWD`
set zipfile = $PWD/$dirname.zip
if ( -e $zipfile) then
  echo "Renaming old $dirname.zip file to $dirname.zip~"
  mv -f $zipfile $zipfile~
endif

cd ..

# zip options:
# -o make zipfile as old as latest entry
# -r recurse into directories
# -x '...' exclude files
if ( "$1" == "zipall" ) then
  zip -or $zipfile $dirname
else 
  zip -or $zipfile $dirname \
    -x '*~' -x '*.mod' -x '*.exe' -x '*.o' -x '*.log' \
    -x '*/Makefile.m' -x '*/workdir*/?*' -x '*.inc_old'
endif

# Symbolic links are now included as the whole files in the zip archive.
# However, some links are internal links (i.e. between directories that
# are both included in the zip file). They must be stored as links. This
# is done by overwriting them in the zip file:

# zip options:
# -y store symbolic links as the link instead of the referenced file
find $dirname/jvpp/dat_m17/spectra   -type l | xargs zip -oy $zipfile
find $dirname/jvpp/dat_lit/spectra   -type l | xargs zip -oy $zipfile
find $dirname/jvpp/dat_m17/hardcoded -type l | xargs zip -oy $zipfile
find $dirname/jvpp/dat_lit/hardcoded -type l | xargs zip -oy $zipfile
zip -oy $zipfile $dirname/messy_cmn_photol_mem.f90
zip -oy $zipfile $dirname/messy_jval_jvpp.inc
zip -oy $zipfile $dirname/mo_netcdf.f90

if ( "$2" == "verbose" ) then
  echo;echo "The zipfile has been created:"
  ls -la $zipfile
endif

exit
