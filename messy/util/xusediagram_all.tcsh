#!/bin/tcsh -f
    
#  op_pj_20180314+
set smlist = (`grep -i LOGICAL messy/smcl/messy_main_switch.f90 | grep USE | awk '{print tolower($3)}' | sort | uniq | sed 's|use_||g'`)

foreach sm ($smlist)
  ./messy/util/xusediagram.tcsh $sm
end

cd workdir
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=ALL.pdf *_use.pdf
cd ..

exit 0
#  op_pj_20180314-

# ./messy/util/xusediagram.tcsh a2o
# ./messy/util/xusediagram.tcsh aeropt
# ./messy/util/xusediagram.tcsh airsea
# ./messy/util/xusediagram.tcsh airtraf
# ./messy/util/xusediagram.tcsh bioburn
# ./messy/util/xusediagram.tcsh bufly
# ./messy/util/xusediagram.tcsh ch4
# ./messy/util/xusediagram.tcsh chemglue
# ./messy/util/xusediagram.tcsh cloud
# ./messy/util/xusediagram.tcsh cloudj
# ./messy/util/xusediagram.tcsh contrail
# ./messy/util/xusediagram.tcsh convect
# ./messy/util/xusediagram.tcsh cvtrans
# ./messy/util/xusediagram.tcsh d14co
# ./messy/util/xusediagram.tcsh ddep
# ./messy/util/xusediagram.tcsh dissoc
# ./messy/util/xusediagram.tcsh diumod
# ./messy/util/xusediagram.tcsh dradon
# ./messy/util/xusediagram.tcsh drydep
# ./messy/util/xusediagram.tcsh e4chem
# ./messy/util/xusediagram.tcsh e5vdiff
# ./messy/util/xusediagram.tcsh gec
# ./messy/util/xusediagram.tcsh gmxe
# ./messy/util/xusediagram.tcsh gwave
# ./messy/util/xusediagram.tcsh h2o
# ./messy/util/xusediagram.tcsh h2oiso
# ./messy/util/xusediagram.tcsh hd
# ./messy/util/xusediagram.tcsh import
# ./messy/util/xusediagram.tcsh jval
# ./messy/util/xusediagram.tcsh jvst
# ./messy/util/xusediagram.tcsh lnox
# ./messy/util/xusediagram.tcsh m7
# ./messy/util/xusediagram.tcsh made
# ./messy/util/xusediagram.tcsh made3
# ./messy/util/xusediagram.tcsh mecca
# ./messy/util/xusediagram.tcsh megan
# ./messy/util/xusediagram.tcsh mlocean
# ./messy/util/xusediagram.tcsh mmforce
# ./messy/util/xusediagram.tcsh mpiom
# ./messy/util/xusediagram.tcsh msbm
# ./messy/util/xusediagram.tcsh mtskip
# ./messy/util/xusediagram.tcsh ncregrid
# ./messy/util/xusediagram.tcsh o3orig
# ./messy/util/xusediagram.tcsh offemis
# ./messy/util/xusediagram.tcsh offlem
# ./messy/util/xusediagram.tcsh onemis
# ./messy/util/xusediagram.tcsh onlem
# ./messy/util/xusediagram.tcsh oracle
# ./messy/util/xusediagram.tcsh orogw
# ./messy/util/xusediagram.tcsh photo
# ./messy/util/xusediagram.tcsh plumegas
# ./messy/util/xusediagram.tcsh ptrac
# ./messy/util/xusediagram.tcsh ptracini
# ./messy/util/xusediagram.tcsh qbo
# ./messy/util/xusediagram.tcsh rad
# ./messy/util/xusediagram.tcsh rad4all
# ./messy/util/xusediagram.tcsh radjimt
# ./messy/util/xusediagram.tcsh relax
# ./messy/util/xusediagram.tcsh s4d
# ./messy/util/xusediagram.tcsh satsims
# ./messy/util/xusediagram.tcsh scalc
# ./messy/util/xusediagram.tcsh scav
# ./messy/util/xusediagram.tcsh scout
# ./messy/util/xusediagram.tcsh sedi
# ./messy/util/xusediagram.tcsh sorbit
# ./messy/util/xusediagram.tcsh spacenox
# ./messy/util/xusediagram.tcsh spe
# ./messy/util/xusediagram.tcsh surface
# ./messy/util/xusediagram.tcsh tagging
# ./messy/util/xusediagram.tcsh tbudget
# ./messy/util/xusediagram.tcsh timepos
# ./messy/util/xusediagram.tcsh tnudge
# ./messy/util/xusediagram.tcsh trexp
# ./messy/util/xusediagram.tcsh tropop
# ./messy/util/xusediagram.tcsh trsync
# ./messy/util/xusediagram.tcsh ubcnox
# ./messy/util/xusediagram.tcsh vahr
# ./messy/util/xusediagram.tcsh vaxtra
# ./messy/util/xusediagram.tcsh vertdiff
# ./messy/util/xusediagram.tcsh vertex
# ./messy/util/xusediagram.tcsh viso
#
#cd workdir
#pdftk *_use.pdf cat output ALL.pdf
