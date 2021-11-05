set term postscript color
set grid
set output 'm7_plot.ps'
# gnuplot Version 3.8i and newer use:
# set style line  1 lt  2 lw 1.5 pt  1 ps 1.00
# set style line  2 lt  1 lw 1.0 pt  3 ps 1.25
# set style line  3 lt  4 lw 1.5 pt  5 ps 1.00
# set style line  4 lt  3 lw 1.5 pt  7 ps 1.25
# set style line  5 lt  5 lw 1.5 pt  9 ps 1.00
# set style line  6 lt  6 lw 1.5 pt 11 ps 1.00
# set style line  7 lt  9 lw 1.5 pt 13 ps 1.25
# set style line  8 lt  8 lw 1.5 pt 15 ps 1.00
# set style line  9 lt  7 lw 1.5 pt 17 ps 1.00
# set style line 10 lt 10 lw 1.5 pt 19 ps 1.00
# set style line 11 lt 11 lw 1.5 pt 21 ps 1.00
# set style line 12 lt 12 lw 1.5 pt 22 ps 1.00
# set style line 13 lt -2 lw 1.0 pt  0 ps 1.00

# gnuplot Version 3.7 and older use:
set  line  1 lt  2 lw 1.5 pt  1 ps 1.00
set  line  2 lt  1 lw 1.0 pt  3 ps 1.25
set  line  3 lt  4 lw 1.5 pt  5 ps 1.00
set  line  4 lt  3 lw 1.5 pt  7 ps 1.25
set  line  5 lt  5 lw 1.5 pt  9 ps 1.00
set  line  6 lt  6 lw 1.5 pt 11 ps 1.00
set  line  7 lt  9 lw 1.5 pt 13 ps 1.25
set  line  8 lt  8 lw 1.5 pt 15 ps 1.00
set  line  9 lt  7 lw 1.5 pt 17 ps 1.00
set  line 10 lt 10 lw 1.5 pt 19 ps 1.00
set  line 11 lt 11 lw 1.5 pt 21 ps 1.00
set  line 12 lt 12 lw 1.5 pt 22 ps 1.00
set  line 13 lt -2 lw 1.0 pt  0 ps 1.00

set xdata time
set timefmt "%d.%m.%y%y.%H:%M"
set format x "%d/%m/%y %H/%M"
set format x "%d/%m/%y"
set format x "%H:%M.%d"
set xrange        ["01.07.2001.12:00":"05.07.2001.12:00"]
set xlabel   'Time [01.07.2001.12:00 - 05.07.2001.12:00]'
set xrange =
set xlabel =
set title    "Aerosol Dynamics with M7 (v1.5)" 

#set logscale y 
#unset logscale y 

set ylabel   'T [deg.C]'
plot 'm7_input_1.dat' using 1:2  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:2  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:2  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:2  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:2  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:2  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:2  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:2  ti 'CI'     w linespoints ls 10
set ylabel   'RH [%]'
plot 'm7_input_1.dat' using 1:3  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:3  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:3  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:3  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:3  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:3  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:3  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:3  ti 'CI'     w linespoints ls 10
set ylabel   'P [hPa]'
plot 'm7_input_1.dat' using 1:4  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:4  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:4  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:4  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:4  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:4  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:4  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:4  ti 'CI'     w linespoints ls 10
set ylabel   'd_wet [m]
plot                                                                'm7_1.dat' using 1:14  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:14  ti 'KS'    w linespoints ls 5, 'm7_3.dat' using 1:14  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:14  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:14  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:14  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:14  ti 'CI'     w linespoints ls 10
set logscale y
set ylabel   'r_dry [m]'
plot                                                                'm7_1.dat' using 1:12  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:12  ti 'KS'    w linespoints ls 5, 'm7_3.dat' using 1:12  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:12  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:12  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:12  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:12  ti 'CI'     w linespoints ls 10
set ylabel   'r_wet [m]'
plot                                                                'm7_1.dat' using 1:13  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:13  ti 'KS'    w linespoints ls 5, 'm7_3.dat' using 1:13  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:13  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:13  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:13  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:13  ti 'CI'     w linespoints ls 10
set ylabel   'H2SO4 [molecules kg-1 (air)]'
plot 'm7_input_1.dat' using 1:5  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:5  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:5  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:5  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:5  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:5  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:5  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:5  ti 'CI'     w linespoints ls 10
set ylabel   'SO4 [molecules kg-1 (air)]'
plot 'm7_input_1.dat' using 1:5  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:6  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:6  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:6  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:6  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:6  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:6  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:6  ti 'CI'     w linespoints ls 10
set ylabel   'BC [ug m-3 (air)]'
plot 'm7_input_1.dat' using 1:6  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:7  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:7  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:7  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:7  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:7  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:7  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:7  ti 'CI'     w linespoints ls 10
set ylabel   'OC [ug m-3 (air)]'
plot 'm7_input_1.dat' using 1:7  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:8  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:8  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:8  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:8  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:8  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:8  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:8  ti 'CI'     w linespoints ls 10
set ylabel   'SS [ug m-3 (air)]'
plot 'm7_input_1.dat' using 1:8  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:9  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:9  ti 'KS'     w linespoints ls 5, 'm7_3.dat' using 1:9  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:9  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:9  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:9  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:9  ti 'CI'     w linespoints ls 10
set ylabel   'DU [ug m-3 (air)]'
plot 'm7_input_1.dat' using 1:9  ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:10  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:10  ti 'KS'    w linespoints ls 5, 'm7_3.dat' using 1:10  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:10  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:10  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:10  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:10  ti 'CI'     w linespoints ls 10
set ylabel   'N [1/cm3 (air)]'
plot 'm7_input_1.dat' using 1:10 ti 'Input'     w linespoints ls 13, 'm7_1.dat' using 1:11  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:11  ti 'KS'    w linespoints ls 5, 'm7_3.dat' using 1:11  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:11  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:11  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:11  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:11  ti 'CI'     w linespoints ls 10
#unset logscale y
set logscale y
set ylabel   'H2O [ug m-3 (air)]'
plot                                                                'm7_1.dat' using 1:15  ti 'NS'     w linespoints ls 4, 'm7_2.dat' using 1:15  ti 'KS'    w linespoints ls 5, 'm7_3.dat' using 1:15  ti 'AS'     w linespoints ls 6, 'm7_4.dat' using 1:15  ti 'CS'     w linespoints ls 7, 'm7_5.dat' using 1:15  ti 'KI'     w linespoints ls 8, 'm7_6.dat' using 1:15  ti 'AI'     w linespoints ls 9, 'm7_7.dat' using 1:15  ti 'CI'     w linespoints ls 10

