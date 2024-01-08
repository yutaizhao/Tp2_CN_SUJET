set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'perf_direct.png'

set title 'Latences de résolution de 2 méthodes'
set xlabel 'taille de matrice'
set ylabel 'latence (ns)'

plot 'PERF_SV' with lines title 'dgbsv', \
     'PERF_TRFnTRS' with lines title 'dgbtrf+dgbtrs'