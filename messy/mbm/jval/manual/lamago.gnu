set term postscript eps color
set output "lamago.eps"
set xlabel "solar zenith angle (degree)"
set ylabel "Fcorr"
set dummy sza
#set logscale y

Fcorr1(sza) = exp(19.09*1  *(1-sza/87.5))
Fcorr2(sza) = exp(19.09*1.3*(1-sza/87.5))
Fcorr3(sza) = exp(19.09*1.4*(1-sza/87.5))
Fcorr4(sza) = exp(19.09*1.5*(1-sza/87.5))
Fcorr5(sza) = exp(19.09*2  *(1-sza/87.5))
Fcorr6(sza) = exp(19.09*3  *(1-sza/87.5))
Fcorr7(sza) = exp(19.09*4  *(1-sza/87.5))
Fcorr8(sza) = exp(19.09*5  *(1-sza/87.5))

plot [87.5:93] Fcorr1(sza), Fcorr2(sza), Fcorr3(sza), Fcorr4(sza), Fcorr5(sza), Fcorr6(sza), Fcorr7(sza), Fcorr8(sza)

