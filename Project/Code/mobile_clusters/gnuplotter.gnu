set terminal png size 1200,1000 enhanced font "Helvetica,12" 
set output "fig/0016.png" 
set xrange [0:99]
set yrange [0:99]
set style fill transparent solid 1.0 noborder
unset key
plot "data/0016.txt" u 1:2:3:4 w circles palette
