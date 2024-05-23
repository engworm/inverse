
# set logscale x
set logscale y
set ytics nomirror
set logscale y2

set y2tics 

set format y "%.0E"
set format y2 "%.0E"

set xlabel "Number of Iteration"
set ylabel "Cost"
set y2label "Condition Number"

set key font"Arial,20"

plot "grav/dp/log/data.txt" with lp axis x1y1 title "Cost (Gravity) " pt 6 ps 1 lc "web-blue", \
     "pote/dp/log/data.txt" with lp axis x1y1 title "Cost (Potential) " pt 6 ps 1 lc "red", \

     # "grav/log/data.txt" u 1:3 with lp axis x1y2 title "Condition Number (Gravity)" pt 70 ps 1 lc "web-blue" dt 3, \
     # "pote/log/data.txt" u 1:3 with lp axis x1y2 title "Condition Number (Potential)" pt 70 ps 1 lc "red" dt 3, \

pause -1
