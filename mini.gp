set size square
set xrange [-2:2]
set yrange [-2:2]

set key below
# unset title
# unset key
set tics font "Arial,20"

plot  "grav/dp/log/xm.txt" u 1:2 pt 1 ps 2 lc 1 title "gravity", \
      "pote/dp/log/xm.txt" u 1:2 pt 2 ps 2 title "potential"

# plot  "grav/mp/log/xm.txt" u 1:2 pt 1 ps 2 lc 1 title "gravity(mp)", \
      # "pote/mp/log/xm.txt" u 1:2 pt 2 ps 2 title "potential(mp)"

pause -1
