set term post eps color solid enh
set output "try.eps"
set title "Lorenz Equations"
set xlabel "X"
set ylabel "Y"
plot "result.txt" using 1:2 title "X-Y Phase Space" with lines


