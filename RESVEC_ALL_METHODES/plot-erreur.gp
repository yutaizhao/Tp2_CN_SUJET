set title 'Multiple Data Files Plot'
set xlabel 'X-axis'
set ylabel 'Y-axis'

set xrange [0:200]

plot 'RESVEC_alpha.dat' with line title 'alpha',\
     'RESVEC_jac.dat' with line title 'Jacobi',\
     'RESVEC_gs.dat' with line title 'Gauss-Seidel'

pause -1

