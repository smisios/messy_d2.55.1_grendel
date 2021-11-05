#!/bin/tcsh -f

#set echo verbose

# postprocess e.g. by
# xchecknml.tcsh | grep "without generic" | awk -F ':' '{printf(" - [ ] %s\n",$3)}' | sort | uniq
# or other greps ...

@ rs = 0
@ ipt = 0
@ nc = 0
@ ncc = 0
@ rp = 0
@ nn = 0
@ dnc = 0
@ dncc = 0

set nmllist = (`find messy/nml -type f -print | sort`)
foreach nmlfile (${nmllist})
  set nmlb = `basename $nmlfile`
  if ("$nmlb" == "README") continue
  if ($nmlfile =~ '*xmessy*') then
     echo "ERROR: RUNSCRIPT in namelist directory: $nmlfile"
     @ rs++
     continue
  endif
  if ($nmlfile !~ '*.nml*') then
     echo "ERROR: non-namelist file: $nmlfile"
     @ nn++
     continue
  endif
  #echo $nmlfile
  #sed 's|\!.*$||g' ${nmlfile}
  set rg = `sed -E 's|\!'\.\*\$'||g' ${nmlfile} | grep -i '&regrid'`
  if ("${rg}" != "") then
     set bn = `basename ${nmlfile}`
     #echo $bn
     set iflist = (`sed -E 's|\!'\.\*\$'||g' ${nmlfile} | grep -i infile | awk -F "=" '{print $2}' | sed 's|[\"'\'',]||g' | tr ' ' '\n'`)
     foreach if ($iflist)
	set p1 = `echo $if | awk -F '/' '{print $1}'`
	#echo $p1
	if ($p1 == '$INPUTDIR_CESM1') continue
	if ($p1 == '$INPUTDIR_ICON') continue
	if ($p1 == '$INPUTDIR_ANA') continue
	if ($p1 == '$RESTART_DIR') continue
	if ($p1 != '$INPUTDIR_MESSY') then
	   echo "ERROR: infile without generic path: ${nmlfile}: ${if}"
	   @ ipt++
	endif
     end

     if (${bn} =~ "*tracer*") then
	### no more test required
     else
	### additional tests:
	### - relative path
	### - basename of nml and nc file identical?
	### - naming convention of files
	#echo IMPORT: $nmlfile
	# - basename of nml and nc file
	set bnnm = `basename ${nmlfile} .nml`
	set qdnnm = `dirname ${nmlfile}`
	set dnnm = `echo ${qdnnm} | sed 's|^.*/import/||g'`
	set dp = `echo $qdnnm | awk -F '/' '{printf("%s/%s/%s/%s\n",$1,$2,$3,$4)}'`
	#echo $dnmm
	foreach if ($iflist)
	  #
	  set p1 = `echo $if | awk -F '/' '{print $1}'`
	  if ($p1 == '$INPUTDIR_CESM1') continue
	  if ($p1 == '$INPUTDIR_ICON') continue
	  if ($p1 == '$INPUTDIR_ANA') continue
	  if ($p1 == '$RESTART_DIR') continue
	  #
	  set bnnc = `basename ${if} .nc`
	  set qdnnc = `dirname ${if}`
	  set dnnc = `echo ${qdnnc} | sed 's|$INPUTDIR_MESSY/||g'`
	  if ("${bnnc}" != "${bnnm}") then
	     echo "ERROR: basename naming convention: ${nmlfile}: ${bnnc}"
	     @ nc++
	     if ("$dp" == "messy/nml/DEFAULTS/import") then
		@ dnc++
	     endif
	  endif
	  if ("${dnnc}" != "${dnnm}") then
	     echo "ERROR: relative path mismatch: ${nmlfile}: ${if}"
	     @ rp++
	  endif
	  set nof = `echo ${bnnm} | awk -F '_' '{print NF}'`
	  #set nod = `echo ${bnnm} | awk -F '-' '{print NF}'`
	  #if ( ($nof != 6) || ($nod !=2) ) then
	  if ($nof != 6) then
	     echo "ERROR: file naming convention: ${nmlfile}"
	     @ ncc++
	     if ("$dp" == "messy/nml/DEFAULTS/import") then
		@ dncc++
	     endif
	  endif
	end
     endif
  endif
end

echo '------------------------------------'
echo '<*> SEVERE ERRORS                   '
echo '------------------------------------'
echo "<*> RUNSCRIPTS detected       : $rs"
echo "    NON-NAMELIST FILES        : $nn"
echo "<*> NON-GENERIC PATH IN infile: $ipt"
echo "    NAMING CONVENTION nml=nc  : $nc"
echo "    NAMING CONVENTION nml=nc  : $dnc (DEFAULTS)"
echo "    NAMING CONVENTION         : $ncc"
echo "    NAMING CONVENTION         : $dncc (DEFAULTS)"
echo "    RELATIVE PATH MISMATCH    : $rp"
echo '------------------------------------'

# count non-acceptable errors
@ error = $rs + $ipt

if ($error > 0) then
   exit 1
endif

exit 0
