#!/usr/bin/tcsh -f

# use it in background:
# nohup ./regrid.tcsh < /dev/null >&! /ptmp/bkern/regrid.log &

# set echo verbose

set INPUTDIR = "/work/mm0062/b302040/analy/ao_qctm_ocx"
set OUTPUTDIR = "/work/mm0062/b302040/analy/ao_qctm_ocx/regrid"
set TEMPLATENML = "tracer_om.nml"
#set NCFILE = "tracer_om"
set NCFILE = "_ave"
set A2O = `pwd`

# remove trailing slash (if there is one)

set INPUTDIR = `echo ${INPUTDIR} | sed 's|/$||g'`
set OUTPUTDIR = `echo ${OUTPUTDIR} | sed 's|/$||g'`

# create output dir (if not existing)

mkdir -p ${OUTPUTDIR}

rm -f regrid_tmp.nml

foreach FILE (`/bin/ls ${INPUTDIR}/*${NCFILE}.nc`)

  set FILENAME = `echo ${FILE} | sed 's|.*/||g;'`

  set OUTPUTFILE = ${OUTPUTDIR}/${FILENAME}

  set TIMESTEPS = `ncdump -h $FILE  | grep 'time = UNLIMITED' | sed 's|\(.*\)// (||g;s| .*||g;'`

  cat ${TEMPLATENML} | sed 's|$INFILE|'${FILE}'|g;s|$OUTFILE|'${OUTPUTFILE}'|g;s|$NT|'${TIMESTEPS}'|g;' > regrid_tmp.nml

  ${A2O}/a2o.exe regrid_tmp.nml

  rm -f regrid_tmp.nml

end
