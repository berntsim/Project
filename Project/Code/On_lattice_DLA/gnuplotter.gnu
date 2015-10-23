set terminal png size 600,500 enhanced font "Helvetica,12" 
set output "figures/test.png" 
set xlabel "x-dim" 
set ylabel "y-dim"
set xrange [0:75]
set yrange [0:75]
plot "data/test.txt" with points pointtype 15
