set title 'Historique de convergence'
set xlabel 'Nombre Itérations'
set ylabel 'Norme du résidu relatif'

set xrange [0:200]

plot 'RESVEC_alpha.dat' with line title 'alpha'

pause -1

