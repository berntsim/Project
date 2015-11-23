set terminal png size 1200,1000 enhanced font "Helvetica,12" 
set output "fig/000.100396.png" 
set xrange [0:399]
set yrange [0:399]
set style fill transparent solid 1.0 noborder
unset key
plot "data/000.100396.txt" u 1:2:3:4 w circles palette
