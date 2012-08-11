set term post eps color solid enh
set output "try.eps"
set title "Lorenz Equations"
set xlabel "time step"
set ylabel "Y"
plot "result.txt" using 2 title "Y" with lines


