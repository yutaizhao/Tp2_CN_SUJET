set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'perf_iterative.png'

set title 'Latences de résolution de méthodes Richardson'
set xlabel 'taille de matrice'
set ylabel 'latence (ns)'

plot 'PERF_ALPHA_time' with lines title 'alpha', \
     'PERF_JAC_time' with lines title 'Jacobi', \
     'PERF_GS_time' with lines title 'Gauss-Seidel'