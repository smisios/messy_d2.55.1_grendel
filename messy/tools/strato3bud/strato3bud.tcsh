#! /bin/tcsh -f

set exp_name = M2TS2000_01
set datapath = /work/bd0617/b302019/TS2000/${exp_name}/tr_O3_sbud
set exepath = ../../../bin



foreach year (2000)

#foreach month (`seq -w 1 1 12`)
foreach month (01)

 rm -f tr_O3_sbud.nc
 ln -s ${datapath}/*${year}${month}*.nc tr_O3_sbud.nc

 ${exepath}/strato3bud.exe

 if (-e stratO3bud.nc) then
    mv -f stratO3bud.nc ${exp_name}_stratO3bud_${year}${month}.nc
 else
    echo "ERROR: file stratO3bud.nc was not created"
    exit 1
 endif

end

end

exit 0
