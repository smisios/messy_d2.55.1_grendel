#! /bin/tcsh -f
# Time-stamp: <2019-06-07 14:16:14 sander>
# _zipcaaba.tcsh: create a zip file of caaba code

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

# define links before leaving current directory:
source _internal_links.tcsh

cd ..

# zip options:
# -o make zipfile as old as latest entry
# -r recurse into directories
# -x '...' exclude files
if ( "$1" == "zipall" ) then
  zip -or $zipfile $dirname
else 
  zip -or $zipfile $dirname \
    -x '*~' -x '*.mod' -x '*.exe' -x '*.o' -x '*.a' -x $dirname/'caaba_*.nc' \
    -x '*.log' -x '*.old' -x '*/ferret.jnl' -x '*.zip' -x '*.tar' \
    -x '*.ps' -x '*.pyc' -x '*/__pycache__/*' \
    -x '*.aux' -x '*.bbl' -x '*.toc' -x '*.blg' \
    -x '*/Makefile.m' -x '*_e5.inc' -x '*e4chem*' \
    -x '*/output/?*' -x $dirname/'.git/*' \
    -x '*/tmp/*'
endif

# Symbolic links are now included as the whole files in the zip archive.
# However, some links are internal links (i.e. between directories that
# are both included in the zip file). They must be stored as links. This
# is done by overwriting them in the zip file.
foreach link ($internal_links)
  # zip option -y to store symbolic links as the link instead of the referenced file
  zip -oy $zipfile $dirname/$link
end

if ( "$2" == "verbose" ) then
  echo;echo "The zipfile has been created:"
  ls -la $zipfile
endif

exit
