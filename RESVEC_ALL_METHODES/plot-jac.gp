set title 'Historique de convergence'
set xlabel 'Nombre Itérations'
set ylabel 'Norme du résidu relatif'

set xrange [0:200]

plot 'RESVEC_jac.dat' with line title 'Jacobi'

pause -1

