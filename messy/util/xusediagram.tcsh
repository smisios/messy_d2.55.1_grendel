#!/bin/tcsh -f

##############################################################################
### This script creates a diagram of the Fortran USE - dependencies
### for a specifid MESSy submodel. It utilizes the graphviz software
### (http://www.graphviz.org/) for visualisation.
### The script requires the depend.mk (make dependencies), more specifically
### echam5/__*/depend.mk and messy/smcl/depend.mk (see input below).
### It needs to be called (after configure and gmake) from the MESSy 
### $BASEDIR with a submodel name as parameter.
### The script will operate in the subdirectory workdir to prevent
### results from being zipped.
##############################################################################
### The original script was written by Rolf Sander (MPIC) and has been
### adapted for common usage by Patrick Joeckel (DLR).
### version 1.0 of 2016-08-15
##############################################################################

### for debugging ...
#set echo verbose

##############################################################################

if ($1 == "" ) then
  echo "Usage: `basename $0` <SUBMODEL>"
  exit 1
endif

# check, if dot (graphviz) is available
which dot >&! /dev/null
if ( $? != 0 ) then
   echo "Error: The 'dot' executable of the graphviz software is not in your PATH\!"
   exit 1
endif

set SM = `echo $1 | awk '{print toupper($1)}'`

##############################################################################

set infile = workdir/${SM}_depend.mk
set outfile = workdir/${SM}_use.dot
set pdffile = workdir/${SM}_use.pdf

if (-e $infile) then
   rm -f $infile
endif
touch $infile

if (-e $outfile) then
   rm -f $outfile
endif
touch $outfile

if (-e $pdffile) then
   rm -f $pdffile
endif

##############################################################################

foreach input (echam5/__*/depend.mk messy/smcl/depend.mk)

awk '{for(i=2;i<=NF;i++) {print $1": "$i}}' $input \
  | grep -i _${SM}\[_.\] >> $infile

# mz_rs_20180118+
# add MECCA002, MECCA003 etc. from Polymecca:
if ($SM == 'MECCA') then
  awk '{for(i=2;i<=NF;i++) {print $1": "$i}}' $input \
    | grep -i _${SM}\[0-9\]\[0-9\]\[0-9\]\[_.\] >> $infile
endif
# mz_rs_20180118-
  
end

##############################################################################

echo 'digraph USEDIAGRAM {'                       >> $outfile
echo '  clusterrank=none;'                        >> $outfile
#echo '  clusterrank=local;'                       >> $outfile
#echo '  rankdir=TB;'                              >> $outfile
echo '  rankdir=RL;'                              >> $outfile
echo '  margin = 0.25;'                           >> $outfile
#echo '  page = "183, 217";'                       >> $outfile
#echo '  page = "366, 434";'                       >> $outfile
#echo '  size = "1,2";'                            >> $outfile
#echo '  size = "16.5, 11.7";'                     >> $outfile
#echo '  rotate=90;'                               >> $outfile
#echo '  ratio=0.707;'                             >> $outfile
echo '  ratio="fill";'                            >> $outfile
echo '  graph [fontname = "Helvetica-Bold",'      >> $outfile
echo '    label = "USE-Diagram of '$SM'",'        >> $outfile
echo '    fontsize = 36];'                        >> $outfile
#echo '    fontsize = 36,'                         >> $outfile
#echo '    size = "6,6" ];'                        >> $outfile
echo '  node[shape=box,'                          >> $outfile
echo '    fontname = "Helvetica-Bold",'           >> $outfile
echo '    fontsize=28];'                          >> $outfile

##############################################################################


##############################################################################

cat ${infile} | sed 's| [^ ]*/| |g' | sed 's|:||g' | sed 's| |\n|g' | sed 's|\.[^.]*$||g' | sort | uniq > workdir/tmp_allfiles
#echo "\nworkdir/tmp_allfiles:" ; cat workdir/tmp_allfiles

cat ${infile} | sed 's|\.o:.*||g' | sort | uniq > workdir/tmp_leftfiles
#echo "\nworkdir/tmp_leftfiles:" ; cat workdir/tmp_leftfiles

# BMIL:
grep -iE '_bi$' workdir/tmp_allfiles > workdir/tmp_bmilfiles
grep -iE '_main_control_' workdir/tmp_allfiles >> workdir/tmp_bmilfiles
#echo "\nworkdir/tmp_bmilfiles:" ; cat workdir/tmp_bmilfiles

# SMIL:
grep -iE '_si$' workdir/tmp_allfiles > workdir/tmp_smilfiles
#echo "\nworkdir/tmp_smilfiles:" ; cat workdir/tmp_smilfiles
#grep -viE 'messy_.*_.*_si' workdir/tmp_smilfiles > workdir/tmp_smilfilessamerank
grep -iE 'messy_.*_.*_si' workdir/tmp_smilfiles > workdir/tmp_smilfilessamerank
#echo "\nworkdir/tmp_smilfilessamerank:" ; cat workdir/tmp_smilfilessamerank

# SMILMEM:
grep -iE '_mem_' workdir/tmp_smilfiles > workdir/tmp_smilmemfiles
#echo "\nworkdir/tmp_smilmemfiles:" ; cat workdir/tmp_smilmemfiles

# SMCL:
grep -ivE '_bi$' workdir/tmp_allfiles | grep -ivE '_si$' | grep -ivE '_e5$' | grep -ivE '_main_' > workdir/tmp_smclfiles
#echo "\nworkdir/tmp_smclfiles:" ; cat workdir/tmp_smclfiles

# SMCLMAIN:
#grep -iE '_main_' workdir/tmp_smclfiles > workdir/tmp_smclmainfiles
grep -ivE '_bi$' workdir/tmp_allfiles | grep -ivE '_si$' | grep -ivE '_e5$' | grep -iE '_main_' > workdir/tmp_smclmainfiles
#echo "\nworkdir/tmp_smclmainfiles:" ; cat workdir/tmp_smclmainfiles

# E5:
grep -iE '_e5$' workdir/tmp_allfiles | grep -ivE '_main_' > workdir/tmp_e5files

##############################################################################

echo '\nsubgraph cluster_BMIL {'           >> $outfile
echo 'label = "BMIL (SMIL)"'               >> $outfile
echo 'node[shape=box,'                     >> $outfile
echo '  color=red,'                        >> $outfile
echo '  style="filled,rounded"];'          >> $outfile
sed 's|$|;|g' workdir/tmp_bmilfiles        >> $outfile
echo '}'                                   >> $outfile

# ----------------------------------------------------------------------------

echo '\nsubgraph cluster_SMIL {'           >> $outfile
echo 'label = "SMIL"'                      >> $outfile
echo 'node[shape=box,'                     >> $outfile
echo '  color=yellow,'                     >> $outfile
echo '  style="filled,rounded"];'          >> $outfile
sed 's|$|;|g' workdir/tmp_smilfiles        >> $outfile
echo '}'                                   >> $outfile

echo "{rank=same; "                        >> $outfile
cat workdir/tmp_smilfilessamerank          >> $outfile
echo "}"                                   >> $outfile

# ----------------------------------------------------------------------------

echo '\nsubgraph cluster_SMILMEM {'        >> $outfile
echo 'label = "SMILMEM"'                   >> $outfile
echo 'node[shape=oval,'                    >> $outfile
echo '  color=yellow,'                     >> $outfile
echo '  style="filled"];'                  >> $outfile
sed 's|$|;|g' workdir/tmp_smilmemfiles     >> $outfile
echo '}'                                   >> $outfile

# ----------------------------------------------------------------------------

echo '\nsubgraph cluster_SMCL {'           >> $outfile
echo 'label = "SMCL"'                      >> $outfile
echo 'node[shape=box,'                     >> $outfile
echo '  style="rounded"];'                 >> $outfile
sed 's|$|;|g' workdir/tmp_smclfiles        >> $outfile
echo '}'                                   >> $outfile

# ----------------------------------------------------------------------------

echo '\nsubgraph cluster_SMCLMAIN {'       >> $outfile
echo 'label = "BMIL (SMCL)"'               >> $outfile
echo 'node[shape=box,'                     >> $outfile
echo '  color=orange,'                     >> $outfile
echo '  style="filled,rounded"];'          >> $outfile
sed 's|$|;|g' workdir/tmp_smclmainfiles    >> $outfile
echo '}'                                   >> $outfile

# ----------------------------------------------------------------------------

echo '\nsubgraph cluster_E5 {'             >> $outfile
echo 'label = "E5"'                        >> $outfile
echo 'node[shape=box,'                     >> $outfile
echo '  color=purple,'                     >> $outfile
echo '  style="filled,rounded"];'          >> $outfile
sed 's|$|;|g' workdir/tmp_e5files          >> $outfile
echo '}'                                   >> $outfile

##############################################################################

# dependencies:

grep -v '^#' ${infile} | \
sed 's|\.mod||g' | \
sed 's|\.o||g' | \
sed 's|\.inc||g' | \
sed 's|:||g' | \
awk '{for (i=2;i<=NF;i++) {s1=substr($i,length($i)-3,4); \
                           if (s1 != ".mod") \
            {n=split($i,a,"/"); printf("\"%s\" -> \"%s\";\n",a[n],$1)} } }' >> $outfile 
echo '}' >> $outfile

##############################################################################

### CLEAN
rm -f workdir/tmp_allfiles
rm -f workdir/tmp_bmilfiles
rm -r workdir/tmp_leftfiles
rm -f workdir/tmp_smclfiles
rm -f workdir/tmp_smclmainfiles
rm -f workdir/tmp_smilfiles
rm -f workdir/tmp_smilfilessamerank
rm -f workdir/tmp_smilmemfiles
rm -f workdir/tmp_e5files 

dot -Tpdf -o $pdffile $outfile

echo ""
echo "To see the result, type:\n"
echo "evince $pdffile &\n"
echo ""

exit
# | grep -v messy_main

##############################################################################
