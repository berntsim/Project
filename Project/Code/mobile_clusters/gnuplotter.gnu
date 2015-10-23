set terminal png size 1200,1000 enhanced font "Helvetica,12" 
set output "fig/1.png" 
set xrange [0:4]
set yrange [0:4]
set style fill transparent solid 1.0 noborder
plot "data/config.txt" with circles fc rgb "red",                    "data/config1.txt" with circles fc rgb "navy" 
