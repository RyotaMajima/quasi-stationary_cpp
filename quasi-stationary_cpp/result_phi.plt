unset multiplot; reset

set xran [-5:5]; set yran [-2:2]
set grid
set xla "x"
set yla "|{/Symbol F}(x)|^2"
se sty fill trans solid 0.5 nobor
pl "./output/output_phi.txt" us 1:2 w l ti "V(x)" lc rgb "black" lw 2, \
 "" us 1:3 w filledcur noti lc rgb "red", \
 "" us 1:4 w l noti lc rgb "red", \
 "" us 1:5 w filledcur noti lc rgb "blue", \
 "" us 1:6 w l noti lc rgb "blue"
