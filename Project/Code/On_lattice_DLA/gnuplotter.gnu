set terminal png size 600,500 enhanced font "Helvetica,12" 
set output "figures/test.png" 
set xlabel "x-dim" 
set ylabel "y-dim"
set xrange [0:51]
set yrange [0:51]
plot "data/test.txt" with points pointtype 16
