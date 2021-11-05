#!/bin/tcsh -f

################################################################
### THIS SCRIPT IS TO UPDATE RESTART DIRECTORIES CREATED WITH
### VERSION <= 2.53.0.1 FOR USAGE WITH VERSION >= 2.53.0.1
###   - input: - link targets of restart_*.nc
###   - ouput: - modified link targets of restart_*.nc
################################################################
### requirements: nco (ncdump, ncatted)
################################################################
### (C) Patrick Joeckel, DLR, May 2017
################################################################
#set echo verbose

### EXAMPLE:
### OLD:
#        :timestep = 720 ;
### NEW:
#        :timestep = 720. ;

### (A) FIND OUT LINK TARGET AND SET FILE NAMES

set target = `ls -la restart_g3b.nc | awk '{print $NF}'`
set tdir   = `dirname $target`
set tfile  = `basename $target`
set tcycle = `echo $tfile | awk -F '_' '{print $2}'`

#echo $target
echo "RESTART: "$tdir
#echo $tfile
echo "CYCLE:   "$tcycle

set here = `pwd`
cd $tdir

foreach file (restart_${tcycle}_*.nc)
  echo $file

  # timestep
  set timestep = `ncdump -h $file | grep :timestep | awk '{print $3}'`
  set lastchar = `echo $timestep | awk '{print l=length($1); substr($0,$l)}'`
  if ("$lastchar" != ".") then
     echo ncatted -h -a timestep,global,o,d,${timestep} $file
          ncatted -h -a timestep,global,o,d,${timestep} $file
  endif

  # #start_date_time
  # set sdt_d = `ncdump -h $file | grep :start_date_time | awk '{print $3}' | sed 's|"||g'`
  # set sdt_t = `ncdump -h $file | grep :start_date_time | awk '{print $4}' | sed 's|"||g'`

  # echo $sdt_d" "$sdt_t

  # set ch = `echo $sdt_t | awk '{print substr($1,7,1)}'`
  # if ("$ch" != ".") then
  #    echo ncatted -h -a start_date_time,global,o,c,\"${sdt_d} ${sdt_t}.000\" $file
  #         ncatted -h -a start_date_time,global,o,c,"${sdt_d} ${sdt_t}.000" $file
  # endif

  # # restart_date_time
  # set sdt_d = `ncdump -h $file | grep :restart_date_time | awk '{print $3}' | sed 's|"||g'`
  # set sdt_t = `ncdump -h $file | grep :restart_date_time | awk '{print $4}' | sed 's|"||g'`

  # echo $sdt_d" "$sdt_t

  # set ch = `echo $sdt_t | awk '{print substr($1,7,1)}'`
  # if ("$ch" != ".") then
  #    echo ncatted -h -a restart_date_time,global,o,c,\"${sdt_d} ${sdt_t}.000\" $file
  #         ncatted -h -a restart_date_time,global,o,c,"${sdt_d} ${sdt_t}.000" $file
  # endif

end

cd $here

exit 0
