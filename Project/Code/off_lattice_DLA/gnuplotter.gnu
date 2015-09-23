set terminal png size 1200,1000 enhanced font "Helvetica,12" 
set output "figures/dla.png" 
set xlabel "x-dim" 
set ylabel "y-dim"
set xrange [0:13000]
set yrange [0:13000]
plot "data/dla.txt" with circles fc rgb "navy" fill solid 
