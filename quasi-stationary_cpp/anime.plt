unset multiplot; reset

TN = 500
set grid
set xran [-5:10]
set yran [-2:3]
set xla "x"
set yla "|{/Symbol y}(x,t)|^2"

do for[i=0:TN]{
 plot sprintf("./output/timeEvo/output%03d.txt", i) us 1:2 w l noti lc rgb "red" lw 2,\
  "" us 1:3 w l tit "V(x)" lc rgb "blue" lw 2
 pause 0.0005
}
