set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'Hist_conv_gs.png'

set title 'Historique de convergence'
set xlabel 'Nombre Itérations'
set ylabel 'Norme du résidu relatif'

set xrange [0:150]

plot 'RESVEC_gs.dat' with line title 'Gauss-Seidel'


