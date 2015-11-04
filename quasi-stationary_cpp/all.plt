unset multiplot; reset
TN = 500
L = 40
t_end = 100.0
set xran [-5:10]
set yran [-2:3]
set grid
set xla "x"
set yla "|{/Symbol y}(x,t)|^2"
plot for[i=0:TN:50] sprintf("./output/timeEvo/output%03d.txt", i) us 1:2 w l \
 ti sprintf("t = %.2f", t_end*i/TN),\
  "" us 1:3 tit "V(x)" w l lc rgb "blue" lw 2
