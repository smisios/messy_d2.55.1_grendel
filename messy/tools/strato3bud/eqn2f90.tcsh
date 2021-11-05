#!/bin/tcsh -f

#set echo verbose

set infile = ../../mbm/caaba/mecca/mecca.eqn

# note: no spaces in summation of several equations ...

set list = (\
k_NO_HO2    G3201 \
k_NO_CH3O2  G4104 \
k_O3_O3P    G1003 \
k_O1D_H2O   G2111 \
k_NO2_O3P   G3105 \
k_O3_HO2    G2107 \
k_OH_O3     G2104 \
k_HO2_O3P   G2106 \
k_O3_H      G2101 \
k_OH_O3P    G2103 \
k_ClO_O3P   G6101 \
k_ClO_HO2   G6204 \
k_BrO_O3P   G7101 \
k_BrO_ClO   G7603b+G7603c \
)
# qqq: G7603a not used for k_BrO_ClO !?

################################################################
# special case:
# k_ClO_ClO    = k_3rd_iupac(temp,cair,2.E-32,4.,1.E-11,0.,0.45) ! <G6102d>

# get line number
@ lno = `awk '{if (toupper($1) == "K_CLO_CLO") print NR}' $infile`
#
# check (and eliminate) continuation lines (with & at end)
set final = 0
while ($final == 0)
  set line = `awk '{if (NR == '$lno') print}' $infile`
  echo -n $line | sed 's|temp|temp(ilo,ila,ilev,it)|g' | sed 's|&||g'
  set final = `echo $line | awk '{if (index($0,"&") > 0) {print 0} else {print 1}}'`
  @ lno ++
end
echo
################################################################

set n = ${#list}
#echo $n

@ cnt = 0
while ($cnt < $n)
  @ cnt ++
  set var = $list[$cnt]
  @ cnt ++
  set eqn = $list[$cnt]

  set eqnlist = (`echo $eqn | sed 's|+| |g'`)
  @ ec = 0
  foreach eqn ($eqnlist)
      @ ec ++

      # echo $var $eqn
      #grep ${eqn} $infile | awk -F ':' '{print $2}' | sed 's|;.*$||g' | sed 's|{.*}||g' | sed 's|temp|temp(ilo,ila,ilev,it)|g'

      set rc = "` grep \<${eqn}\> $infile | awk -F ':' '{print "\$2"}' | sed 's|;.*"\$"||g' | sed 's|{.*}||g' | sed 's|temp|temp(ilo,ila,ilev,it)|g' | sed 's| ||g'`"

      if ($ec == 1) then

         if ("$rc" == "") then
            echo "$var = 0.0 ! EQUATION $eqn NOT AVAILABLE"
         else
            echo "$var = $rc ! $eqn"
         endif

      else

         if ("$rc" == "") then
            echo "$var = $var + 0.0 ! EQUATION $eqn NOT AVAILABLE"
         else
            echo "$var = $var + $rc ! + $eqn"
         endif

      endif
   end

end

exit 0


