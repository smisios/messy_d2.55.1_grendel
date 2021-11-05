#! /bin/sh
case $GRID in
  TP01)
    #formally a 0.1deg grid
    reflon='1.45'
    reflat='55'
    np=896
    ie=$(( ${np} * 4 + 2  ))
    je=1800
    ke=${ke-80}
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  TP6M) # how to account for weghting with cos of latitude?
    #formally a 6min grid
    reflon='1.45'
    reflat='55'
    np=900
    ie=$(( ${np} * 4 + 2  ))
    je=2394
    ke=${ke-80}
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  TP01l3)
    reflon='1.45'
    reflat='55'
    np=896
    ie=$(( ${np} * 4 + 2  ))
    je=1800
    ke=3
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  TP10l22)
    # formally a 1.0 deg grid
    reflon='1.45'
    reflat='55'
    np=90
    ie=$(( ${np} * 4 + 2  ))
    je=192
    ke=22
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  TP04)
    # formally a 0.45deg grid
    reflon='1.45'
    reflat='55'
    np=200
    ie=$(( ${np} * 4 + 2  ))
    je=404
    ke=${ke-40}
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  TP10)
    # formally a 1.0 deg grid
    reflon='1.45'
    reflat='55'
    np=90
    ie=$(( ${np} * 4 + 2  ))
    je=192
    ke=${ke-40}
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  TP40)
    # formally a 4.0 deg grid
    reflon='1.45'
    reflat='55'
    np=22
    ie=$(( ${np} * 4 + 2  ))
    je=48
    ke=${ke-20}
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  TOY)
    # formally a 4.0 deg grid
    reflon='1.45'
    reflat='55'
    np=16
    ie=$(( ${np} * 4 + 2  ))
    je=36
    ke=${ke-10}
    tp=true
    periodic_bound='.true.'
    FRMT="F6.0"
    ;;
  GR60)
    ie=60
    je=50
    ke=${ke-20}
    tp=false
    rlat1=72
    rlon1=-40
    rlat2=-84
    rlon2=60
    phi=1.49
    periodic_bound='.true.'
    FRMT="F5.0"
    ;;
  GR30)
    tp=false
    ie=122
    je=101
    ke=${ke-40}
    tp=false
    periodic_bound='.true.'
    echo "FIXME rlat1 rlon1 rlat2 rlon2 phi are no longer known for grid $GRID" >&2
    FRMT="F5.0"
    ;;
  GR15)
    tp=false
    ie=256
    je=220
    ke=${ke-40}
    tp=false
    periodic_bound='.true.'
    echo "FIXME rlat1 rlon1 rlat2 rlon2 phi are no longer known for grid $GRID" >&2
    FRMT="F5.0"
    ;;
  # FIXME: the following incomplete grid configurations are currently
  # only suitable for plotting
  NASH15)
    ie=241
    je=164
    ;;
  T31)
    ie=96
    je=48
    ;;
  1DEG)
    ie=360
    je=180
    ;;
  *)
    echo "Unknown grid configuration ${GRID}." >&2
esac

case $ke in
  3)
    depths="20. 20. 5000."
    ;;
  10)
    depths="40.  50.  90.  160. 270. 390. 550. 900. 1300. 2300."
    ;;
  20)
    depths="20.  20.  20.  30.  40.  50.  70.  90. 120.  150.
           180. 210. 250. 300. 400. 500. 600. 700. 900. 1400."
    ;;
  22)
    depths="50.  50.  78.  54. 110.  64. 148.  78. 196. 102.
           250. 138. 320. 186. 406. 254. 514. 346. 648. 474.
           816. 642."
    ;;
  40)
    depths="12.  10.  10.  10.  10.  10.  13.  15.  20.  25.
            30.  35.  40.  45.  50.  55.  60.  70.  80.  90.
           100. 110. 120. 130. 140. 150. 170. 180. 190. 200.
           220. 250. 270. 300. 350. 400. 450. 500. 500. 600."
    ;;
  80)
    depths="12.  10.  10.  10.  10.  10.  10.  11.  11.  12.
            13.  13.  14.  14.  15.  16.  16.  17.  18.  19.
            20.  21.  21.  22.  24.  25.  26.  27.  28.  29.
            31.  32.  34.  35.  37.  39.  40.  42.  44.  46.
            48.  50.  53.  55.  58.  60.  63.  66.  69.  72.
            76.  79.  83.  87.  91.  95.  99. 104. 108. 113.
           119. 124. 130. 136. 142. 149. 155. 163. 170. 178.
           186. 195. 204. 213. 223. 233. 244. 255. 267. 279."
    ;;
esac
if [ x${depths+set} = xset ]; then
  depths=`echo $depths`
fi

for i in tp ie je ke rlat1 rlon1 rlat2 rlon2 phi \
  reflon reflat np periodic_bound FRMT depths ; do
  if eval test \"x\$\{$i+set\}\" = xset ; then
    export $i
  fi
done
