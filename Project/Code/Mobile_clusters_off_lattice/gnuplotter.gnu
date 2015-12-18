set terminal png size 1200,1000 enhanced font "Helvetica,12" 
set output "fig/003000.png" 
set xrange [0:200]
set yrange [0:200]
set style fill transparent solid 1.0 noborder
unset key
plot "data/003000.txt" u 1:2:3:4 w circles palette
