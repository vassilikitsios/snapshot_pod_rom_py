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
set title "POD temporal modes"

set xlabel "time"
set xrange [*:*]

set ylabel "amplitude"
set yrange [*:*]
set format y "%g"

plot "../../1modes/results/POD.temporal_mode_0001.dat" using 1:2 axes x1y1 title "mode1" with lines lt 2 lc 1 lw 4,\
     "../../1modes/results/POD.temporal_mode_0002.dat" using 1:2 axes x1y1 title "mode2" with lines lt 2 lc 2 lw 4,\
     "../results/ROM.temporal_mode_0001.dat" using 1:2 axes x1y1 title "mode1 ROM" with lines lt 1 lc 1,\
     "../results/ROM.temporal_mode_0002.dat" using 1:2 axes x1y1 title "mode2 ROM" with lines lt 1 lc 2

plot "../../1modes/results/POD.temporal_mode_0003.dat" using 1:2 axes x1y1 title "mode3" with lines lt 2 lc 1 lw 4,\
     "../../1modes/results/POD.temporal_mode_0004.dat" using 1:2 axes x1y1 title "mode4" with lines lt 2 lc 2 lw 4,\
     "../results/ROM.temporal_mode_0003.dat" using 1:2 axes x1y1 title "mode3 ROM" with lines lt 1 lc 1,\
     "../results/ROM.temporal_mode_0004.dat" using 1:2 axes x1y1 title "mode4 ROM" with lines lt 1 lc 2

plot "../../1modes/results/POD.temporal_mode_0005.dat" using 1:2 axes x1y1 title "mode5" with lines lt 2 lc 1 lw 4,\
     "../../1modes/results/POD.temporal_mode_0006.dat" using 1:2 axes x1y1 title "mode6" with lines lt 2 lc 2 lw 4,\
     "../results/ROM.temporal_mode_0005.dat" using 1:2 axes x1y1 title "mode5 ROM" with lines lt 1 lc 1,\
     "../results/ROM.temporal_mode_0006.dat" using 1:2 axes x1y1 title "mode6 ROM" with lines lt 1 lc 2

# -----------------------------------------------------
