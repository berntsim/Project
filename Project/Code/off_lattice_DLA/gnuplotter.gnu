f(x) = a*x+b
fit f(x) 'data/loglog.txt' u 1:2 via a,b
set terminal png size 1200,1000 enhanced font "Helvetica,12" 
set output "figures/loglog(1.7248).png" 
set xlabel "log(N) " 
set ylabel "log(R_g)"
set xrange [0:10]
set yrange [0:15]
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f', a, b)
plot "data/loglog.txt" u 1:2, f(x)
