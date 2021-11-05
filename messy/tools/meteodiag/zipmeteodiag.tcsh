#! /bin/tcsh -f
# Time-stamp: <2015-12-08 18:27:32 joec_pa>
# zipcaaba: create a zip file of meteodiag code

# if the current directory is only a link, the following cd command
# will jump into the real directory:
cd `/bin/pwd`

if ( "$1" == "" ) then
  echo "This script shoud be used via the Makefile, e.g.:"
  echo "  gmake zip     --> archive important files"
  echo "  gmake zipall  --> archive all files"
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
# -y store symbolic links as the link instead of the referenced file
# -x '...' exclude files
if ( "$1" == "zipall" ) then
  zip -or $zipfile $dirname
  else 
  zip -or $zipfile $dirname \
    -x '*~' -x '*.mod' -x '*.exe' -x 'F*.f' -x '*.i90' -x '*.optrpt' \
    -x '*.o' -x '*.a'  \
    -x '*.log' -x '*.old' -x '*/ferret.jnl' -x '*.zip' -x '*.tar' \
    -x '*.ps' -x '*.dat' -x '*/Makefile.m' -x $dirname/'output/?*' \
    -x '*/tmp_*' \
    -x '*.aux' -x '*.bbl' -x '*.toc' -x '*.blg' \
    -x $dirname/'temporaryfile*' \
    -x '*/tmp/*' 
endif

echo "\nThe zipfile has been created:"
ls -la $zipfile

exit
