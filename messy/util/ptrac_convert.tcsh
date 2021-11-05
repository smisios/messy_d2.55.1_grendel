#!/bin/tcsh -f

### This script is to convert pre-d2.52k namelists of the PTRAC submodel
### into the new syntax (d2.52k and above).
### Author: Patrick Joeckel, DLR, August 2016
### Note: This is a quick and (very) dirty solution w.r.t. the coding
###       of this script.
###       YOU NEED TO CHECK THE CONVERTED NAMELIST FILES CAREFULLY!

#set echo verbose

if ($1 == '') then
   echo "Usage: "`basename $0`" <ptrac-namelist file>"
   exit 1
endif
set infile = $1


set tmpfile = tmpfile-`date +"%y%m%d%H%M%S"`

cat $infile \
    | sed 's|\!.*||g' \
    | sed 's|( *\([0-9]*\) *)|(\1)|g' \
    | sed 's|( *\([0-9]*\) *,|(\1,|g' \
    | grep -Ev '^ *$' >> $tmpfile

#cat $tmpfile
#rm -fr $tmpfile
#exit

set no_list = ( \
`grep -i C_NAME $tmpfile | awk '{print $1}' | sed 's|(| |g' | sed 's|)| |g' | awk '{print $2}'`)

cat <<EOF
! -*- f90 -*-
!
!##############################################################################
!### SYNTAX FOR BASIC TRACER DEFINITION:
!###
!### TRAC(.) = 'list of tracer sets', 'list of tracer names', unit, \
!###            medium, quantity, type, \
!###            [aerosol model], [mode], [radius], [sigma], [density],
!###
!###     - lists are separated by semicolon
!###     - medium: 
!###                  AIR        = 1 (default)
!###                  AEROSOL    = 2
!###                  CLOUD      = 3
!###                  OCEAN      = 4
!###                  LAKE       = 5
!###                  RIVER      = 6
!###                  LANDICE    = 7
!###                  SEAICE     = 8
!###                  VEGETATION = 9
!###     - quantity:
!###                  AMOUNTFRACTION = 1 ! = MOLAR MIXING RATIO (default)
!###                  NUMBERDENSITY  = 2
!###                  CONCENTRATION  = 3
!###     - type:
!###                  SINGLE  = 0 (default)
!###                  FAMILY  = 1
!###                  ISOTOPE = 2
!###
!###     - for medium == AEROSOL:
!###                  aerosol model = aerosol model for aerosol properties 
!###                  mode          = aersol mode number [1]
!###     - for aerosol model == 'ptrac' (default):
!###                  radius        = mean aerosol radius [m]
!###                  sigma         = sigma of radius distribution [1]
!###                  density       = aerosol density [kg/m^3]
!###
!##############################################################################
!
&CPL
!
EOF

foreach no ($no_list)

set C_SETS = (`grep -i 'C_SETS('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{for(i=1;i<=NF;i++) {if ($i != "") printf("%s;",$i)}; printf("\n")}'`)
  
# echo $C_SETS

 set C_NAME = `grep -i 'C_NAME('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{if ($2 == "") {printf("%s\n",$1)} else {printf("%s_%s\n",$1,$2)}}'`

#echo $C_NAME

 set C_UNIT = `grep -i 'C_UNIT('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{print $1}'`

# echo $C_UNIT

 set I_TYPE = `grep -i 'I_TYPE('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{print $1}'`

# echo $I_TYPE

 set I_MEDIUM = `grep -i 'I_MEDIUM('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{print $1}'`

# echo $I_MEDIUM

 set I_QUANTITY = `grep -i 'I_QUANTITY('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{print $1}'`

#echo $I_QUANTITY

if ($I_MEDIUM == 2) then # AEROSOL

   set CASK_I = (`grep -i 'CASK_I('$no',' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | sed 's|,| |g'`)
   set AEROSOL_MODE = $CASK_I[14]

   set R_AEROSOL_RADIUS = `grep -i 'R_AEROSOL_RADIUS('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{print $1}'`

   set R_AEROSOL_SIGMA = `grep -i 'R_AEROSOL_SIGMA('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{print $1}'`

   set CASK_R = (`grep -i 'CASK_R('$no',' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | sed 's|,| |g'`)
   set R_AERSOL_DENS = $CASK_R[5]

   set CASK_S = (`grep -i 'CASK_S('$no',' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | sed 's|,| |g'`)
   set S_AERSOL_MODEL = $CASK_S[1]

   echo 'TRAC('$no') = '\'$C_SETS\', \'$C_NAME';'\', \'$C_UNIT\', $I_MEDIUM, $I_QUANTITY, $I_TYPE, \'$S_AERSOL_MODEL\', $AEROSOL_MODE, $R_AEROSOL_RADIUS, $R_AEROSOL_SIGMA , $R_AERSOL_DENS,

else

    echo 'TRAC('$no') = '\'$C_SETS\', \'$C_NAME';'\', \'$C_UNIT\', $I_MEDIUM, $I_QUANTITY, $I_TYPE, , , , , ,

endif

end

cat <<EOF
!
!##############################################################################
!### SYNTAX FOR TRACER PROPERTIES DEFINITION:
!###
!### TPROP(.) = 'list of tracer sets', 'list of tracer names', \
!###            'container name', 'container contents',
!###
!###     - lists are separated by semicolon
!###     - container names are case sensitive
!###     - container content is case sensitive
!###     - ON = 1; OFF = 0
!##############################################################################
!
EOF

@ mo = 1

foreach no ($no_list)

 set C_SETS = `grep -i 'C_SETS('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{for(i=1;i<=NF;i++) {if ($i != "") printf("%s;",$i)}; printf("\n")}'`
  
# echo $C_SETS

 set C_NAME = `grep -i 'C_NAME('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{if ($2 == "") {printf("%s\n",$1)} else {printf("%s_%s\n",$1,$2)}}'`

# echo $C_NAME

 set I_MEDIUM = `grep -i 'C_MEDIUM('$no')' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | awk -F ',' '{for(i=1;i<=NF;i++) {if ($i != "") printf("%s;",$i)}; printf("\n")}'`

 set CASK_I = (`grep -i 'CASK_I('$no',' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | sed 's|,| |g'`)

 set cask_i_strings =  ( \
       'advect        '  'convect       '  'vdiff         '  \
       'wetdep        '  'drydep        '  'sedi          '  \
       'scav          '  'mix           '  'force_col     '  \
       'integrate     '  'timefilter    '  'force_init    '  \
       'aerosol_method'  'aerosol_mode  '  'aerosol_sol   '  \
       'aerosol_hetice'  'hori_diff     '  'relaxation    '  \
       'mmd_init      '   \
       'tag_reg_idt   '  'tag_specific  '                    \
       'lateral_bounds'  'initial type  '  'damping       '  \
       'grib table num'  'grib param num'  'charge        '  \
       'mtskip        ' )

  set cask_i_default = ( 1  1  1  0  0  0  0  1  0  1  1  0 \
         2  0  1  0  0  1  0  0  0 \
         1  1  1   -1   -1  0 0 )

 @ cnt = 1
 foreach i ($CASK_I)
   if ( ($cnt == 14) || ($cnt == 9)) then
      @ cnt ++
      continue
   endif
   if ($i != $cask_i_default[$cnt]) then
   echo 'TPROP('$mo') = '\'$C_SETS\', \'$C_NAME\', \'$cask_i_strings[$cnt]\', \'$i\', 
  @ mo++
   endif
  @ cnt ++
 end

 set CASK_R = (`grep -i 'CASK_R('$no',' $tmpfile | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's| *||g' | sed 's|,| |g'`)

 set cask_r_strings = ( \
       'molarmass      ' 'henry          ' 'dryreac_sf     ' \
       'vini           ' 'aerosol_density' )
 set cask_r_default = ( 0.0 0.0 0.0 0.0 0.0 )

 @ cnt = 1
 foreach i ($CASK_R)
   if ($cnt == 5) then
      @ cnt ++
      continue
   endif
   if ($i != $cask_r_default[$cnt]) then
   echo 'TPROP('$mo') = '\'$C_SETS\', \'$C_NAME\', \'$cask_r_strings[$cnt]\', \'$i\', | sed 's|henry|pss|g'
  @ mo++
  endif
  @ cnt ++
 end

end

cat <<EOF
!
/
EOF

rm -fr $tmpfile

exit 0 
