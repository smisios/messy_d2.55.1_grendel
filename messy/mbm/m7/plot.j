#! /bin/tcsh -f
################################################################################
set SYSTEM = `uname`
set editor = "emacs"
set launch = "ls -l"
set gnuplot = gnuplot
set ps2pdf  = ps2pdf
set launch  = gv
set date   = `date +"20%y.%m.%d_%Hhr%M"`
################################################################################
if ($SYSTEM == "Darwin") then
set editor  = "bbedit"
set launch  = open
set ps2pdf  = ps2pdf
set gnuplot = gnuplot
endif
################################################################################
cd output
################################################################################
#$editor    ../m7.plot
set last_data_date = `cat last_data_date`
set first_calculation_date = `cat first_calculation_date`
set last_calculation_date  = `cat last_calculation_date`
if ($last_data_date != $last_calculation_date) then
echo "last data date ($last_data_date) differs from last calculation date ($last_calculation_date)!"
echo "=> the last data date will be used for ploting to make use of the input data for reference."
set last_calculation_date = $last_data_date
endif
cat   ../m7.plot | sed "s|set xrange =|set xrange        ['$first_calculation_date':'$last_calculation_date']|g" > tmp
cat      tmp     | sed "s|set xlabel =|set xlabel   'Time [$first_calculation_date - $last_calculation_date]'|g" > m7.plot
$gnuplot <    m7.plot
$ps2pdf       m7_plot.ps
rm            m7_plot.ps tmp
mv m7_plot.pdf ./m7_v1.5_plot_$date.pdf
$launch        ./m7_v1.5_plot_$date.pdf
################################################################################
cd ..
################################################################################
exit
################################################################################
