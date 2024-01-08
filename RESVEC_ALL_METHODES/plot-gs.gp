set title 'Historique de convergence'
set xlabel 'Nombre Itérations'
set ylabel 'Norme du résidu relatif'

set xrange [0:200]

plot 'RESVEC_gs.dat' with line title 'Gauss-Seidel'

pause -1

