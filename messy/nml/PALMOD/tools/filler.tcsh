#! /bin/tcsh

#set cparts = ( cma cmd cda cdd ima imd ) #ida idd )
 set cparts = ( c c0 m )

 set afile = $1
#set afile = gas_ohvar.eqn
#set afile = gas.tex

rm $0.tmp
foreach cp ( ${cparts} )

  grep '##' ${afile} | sed -e "s/##/${cp}/g" >> $0.tmp

end

cat $0.tmp | sort

rm $0.tmp
