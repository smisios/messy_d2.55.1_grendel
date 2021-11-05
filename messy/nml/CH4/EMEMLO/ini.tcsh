#!/bin/tcsh -f

foreach file ( restart_0001_tracer_gp_m1.nc restart_0001_tracer_gp.nc)

 cp ${file} ${file}-orig
 ncks -O -x -v SF6,AOA,SF6_AOA,SF6_AOAc,SF6_CCMI ${file}-orig ${file}
 ncks -A -v SF6,AOA,SF6_AOA,SF6_AOAc,SF6_CCMI ../../../R1EME-03-oc/save/0024/${file} ${file}

 end

 # set file = restart_0001_tracer_pdef_gp.nc

 # cp ${file} ${file}-orig
 # ncks -O -x -v MP_SF6_x1,MP_AOA_x1,MP_SF6_AOA_x1,MP_SF6_AOAc_x1,MP_SF6_CCMI_x1,MN_SF6_x1,MN_AOA_x1,MN_SF6_AOA_x1,MN_SF6_AOAc_x1,MN_SF6_CCMI_x1 ${file}-orig ${file}

 # ncks -A -v MP_SF6_x1,MP_AOA_x1,MP_SF6_AOA_x1,MP_SF6_AOAc_x1,MP_SF6_CCMI_x1,MN_SF6_x1,MN_AOA_x1,MN_SF6_AOA_x1,MN_SF6_AOAc_x1,MN_SF6_CCMI_x1 ../../../R1EME-03-oc/save/0024/${file} ${file}
 
exit 0
