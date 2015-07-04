# PLOT DATA USING GNUPLOT
#
#> gnuplot plot_output.gp
#
#-------------------------------- Ouput settings ------------------------------
set terminal postscript enhanced color font "Times-Roman" 20 dashlength 3
#set terminal postscript enhanced eps font "Times-Roman" 32 dashlength 3

# ------------------------------------------------------
# functions
phase(x) = x<0 ? x+3.147 : x

# ------------------------------------------------------
set output "output.eps"
set key top right
set key spacing 1.5

# -----------------------------------------------------
set title "POD eignevalues (energies)"
set xlabel "mode number"
set logscale x
set xrange [*:*]
set ylabel "percentage energy"
set logscale y
set yrange [*:*]
set format y "10^{%L}"
plot "../results/POD.eigenvalues.dat" using 1:4 axes x1y1 title "" with linespoints lt 1 lc 1 pt 3

unset logscale y
set format y "%g"

# -----------------------------------------------------
