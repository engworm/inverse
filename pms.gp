# set term lua tikz
# set term tikz size 5cm,5cm
# set output "paper/fig/tg_mpR4.tex"

set pm3d map
set size square
set cbr [0:11]
set palette cubehelix start 3 cycles 0.5 saturation 1 negative


set grid linewidth 2

# set title font "Arial,30"

# set multiplot layout 2,3


set xrange [100:500]
set yrange [100:500]

# set xtics ("-2" 0, "-1.5" 50, "-1" 100, "-0.5" 150, "0" 200, "0.5" 250, "1" 300, "1.5" 350, "2" 400)
# set ytics ("-2" 0, "-1.5" 50, "-1" 100, "-0.5" 150, "0" 200, "0.5" 250, "1" 300, "1.5" 350, "2" 400)

set xtics ("-3" 0, "-2.5" 50, "-2" 100, "-1.5" 150, "-1" 200, "-0.5" 250, "0" 300, "0.5" 350, "1" 400, "1.5" 450, "2" 500, "2.5" 550, "3" 600)
set ytics ("-3" 0, "-2.5" 50, "-2" 100, "-1.5" 150, "-1" 200, "-0.5" 250, "0" 300, "0.5" 350, "1" 400, "1.5" 450, "2" 500, "2.5" 550, "3" 600)

splot "grav/PMS/scatter30.txt" matrix
# splot "pote/PMS/scatter33780.txt" matrix

# splot "PMS/answer/triangle.txt" matrix
# splot "PMS/answer/elliptic.txt" matrix

pause -1
