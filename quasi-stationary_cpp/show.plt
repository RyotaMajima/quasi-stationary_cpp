unset multiplot; reset
set ter tikz size 257mm,364mm standalone
set output "./graph/result.tex"
set multiplot layout 4,1
set grid lw 2

E0 = -1.02
E1 = -0.1575
E2 = 0.5325

set label 1 "$E_{0} = -1.02$" at first E0-0.1, 0.45 font ",14"
set label 2 "$E_{1} = -0.1575$" at fir E1-0.1, 0.45 font ",14"
set label 3 "$E_{2} = 0.5325$" at fir E2-0.1, 0.45 font ",14"
set arrow 1 from first E0,0.4 to E0,0.2 lw 2
set arrow 2 from first E1,0.4 to E1,0.2 lw 2
set arrow 3 from first E2,0.4 to E2,0.05 lw 2
set yran [0:1]
set xla "$\\varepsilon$"
set yla "$|\\Phi_{T}(\\varepsilon)|^{2}$"
set title "peak of energy eigenvalue"
pl "./output/energy.txt" w l noti lc rgb "red" lw 3

unse arrow; unse label
se sty fill trans solid 0.5 nobor

set xran [-5:5]; set yran [-2:2]
set xla "$x$"; set yla "$|\\Phi_{T}(x, \\varepsilon)|^{2}$"

set title "ground state"
pl "./output/output_phi.txt" us 1:2 w l noti lc rgb "black" lw 2, \
 "" us 1:3 w filledcur ti "calc" lc rgb "red", \
 "" us 1:4 w l dt 2 lc rgb "red" noti

se title "first excited state"
pl "" us 1:2 w l noti lc rgb "black" lw 2, \
 "" us 1:5 w filledcur ti "calc" lc rgb "blue", \
 "" us 1:6 w l dt 2 lc rgb "blue" noti

se title "second excited state"
pl "" us 1:2 w l noti lc rgb "black" lw 2, \
 "" us 1:7 w filledcur ti "calc" lc rgb "green", \
 "" us 1:8 w l dt 2 lc rgb "green" noti

unset multiplot
set output
set ter pop