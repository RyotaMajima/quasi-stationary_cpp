unset multiplot; reset

set yran [0:1]
set grid
set xla "{/Symbol e}"
set yla "|{/Symbol f}({/Symbol e})|^2"
pl "./output/energy.txt" us 1:(1*$2) w l noti lc rgb "blue" lw 2
