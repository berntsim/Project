set terminal png size 1200,1200 enhanced font "Helvetica,12" 
set output "fig/0-1.png" 
set xrange [0:10]
set yrange [0:10]
set style fill transparent solid 0.5 noborder
plot "data/before.txt" with circles fc rgb "red",                    "data/after.txt" with circles fc rgb "navy" 
