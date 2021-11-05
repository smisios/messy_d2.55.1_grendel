#! /bin/tcsh -f

# called from xmecca
#
# Purpose:
# postprocess files such as "messy_mecca_kpp.f90" and "messy_mtchem_kpp.f90" 
# and check for lines such as 
#    IF (ind_xyz>0) R02 = RO2 +  C(ind_xyz)
# and
#    INTEGER, PARAMETER,PUBLIC :: ind_xyz = 0
# to remove lines that might contain an "IF ... = ... + C(0)" and
# hence are in conflict with NAG compiler.
#
# Notes:
# - Side effect: gain performance by avoiding unnecessary IFs.
# - RO2 will still be calculated with all indices C(i) with
#   "i" not equal zero.
#
# Known issues:
# 1) The script currently does only work after kp4, since only then
#    all required lines are in messy_<submodel>_kpp.f90 file; without
#    kp4, the lines with
#     - IF conditions are in messy_<submodel>_kpp_rates.f90
#     - ind_ definitions are in messy_<submodel>_kpp_parameters.f90,
#       howeve, without the "PUBLIC" attribute.
# 2) The script is rather hardwired w.r.t. syntax, e.g. upper/lower-case,
#    white spaces, etc. It might be an option to replace it by a more
#    robust python3 script soon.
# 3) It is unclear, if the script works correctly for polyMECCA
# 4) It is unclear, if the script works with the vectorized solver
#    rosenbrock_vec. This implies the question, whether RO2 is required
#    to be vectorized as well.
#

set  inkppf90 = $1

### 1) get species names from corresponding lines

set slist = ( \
`grep -i '^ *IF *( *ind_.*> *0 *)' $inkppf90 \
  | awk -F "[(>]" '{print $2}' \
  | awk -F '_' '{print $2}' \
  | sed 's| ||g'` \
)

### 2) look for each species what value is set for ind_xyz

set no = ()

foreach spec ($slist)

 set n = `\
 grep -i "INTEGER, PARAMETER,PUBLIC :: ind_$spec " $inkppf90 \
  | awk -F " = " '{print $2}'`

 set no = ($no $n)

end

echo $no
echo ${#slist}
echo ${#no}

### 3) an intermediate check
if (${#slist} != ${#no}) then
   echo "ERROR in "$0": number of species found for >>(if ind_xyz > 0)<< and >>(PARAMETER :: ind_xyz)<< mismatch:"
   echo "               ${#slist} != ${#no}"
   exit 1
endif

### 4) copy original file and modify with sed -i ...
# i.e. comment out lines for ind_xyz equals to "0"

cp -f $inkppf90 ${inkppf90}.bkp.old

@ cnt = 1
while ($cnt <= ${#no})
 if (${no[$cnt]} == 0) then
    echo $cnt" "${slist[$cnt]}" <remove>"
    sed -i '/^ *IF *( *ind_'${slist[$cnt]}'.*> *0 *)/d' $inkppf90
    sed -i '/^ *IF *(( *ind_'${slist[$cnt]}'.*> *0 *)/d' $inkppf90
 else
    echo $cnt" "${slist[$cnt]}" <keep>"
    sed -i '/^ *IF *( *ind_'${slist[$cnt]}'.*> *0 *)/s/^.*RO2 =/  RO2 =/' $inkppf90
 endif
@ cnt++
end

############### 2nd part delete statements from F90_INIT ###########
####################################################################
### B 1) get species names from corresponding lines
set s2list = ( \
`grep -i 'IF *((ind_.*> *0*) *' $inkppf90 \
  | awk -F '>' '{print $1}' \
  | awk -F '_' '{print $2}' \
  | sed 's| ||g'` \
)
echo " list2"
echo $s2list

echo ${#s2list}
### B 2) look for each species what value is set for ind_xyz

set no2 = ()

foreach spec ($s2list)

 set n2 = `\
 grep -i "INTEGER, PARAMETER,PUBLIC :: ind_$spec " $inkppf90 \
  | awk -F " = " '{print $2}'`

 set no2 = ($no2 $n2)

end

echo $no2
echo ${#s2list}
echo ${#no2}

### B 3) an intermediate check
if (${#s2list} != ${#no2}) then
   echo "ERROR in "$0": number of species found for >>(if (ind_xyz > 0).AND.<< and >>(PARAMETER :: ind_xyz)<< mismatch:"
   echo "               ${#s2list} != ${#no2}"
   exit 1
endif

### B 4) copy original file and modify with sed -i ...
# i.e. comment out lines for ind_xyz equals to "0"

#cp -f $inkppf90 ${inkppf90}.bkp.old

@ cnt = 1
while ($cnt <= ${#no2})
 if (${no2[$cnt]} == 0) then
    echo $cnt" "${s2list[$cnt]}" <remove>"
    sed -i '/^ *IF *(( *ind_'${s2list[$cnt]}'.*> *0 *)/d' $inkppf90
 else
    echo $cnt" "${s2list[$cnt]}" <keep>"
    #sed -i '/^ *IF *(( *ind_'${s2list[$cnt]}'.*> *0 *)/s/^.*atol/  atol/' $inkppf90
 endif
@ cnt++
end

exit 0
