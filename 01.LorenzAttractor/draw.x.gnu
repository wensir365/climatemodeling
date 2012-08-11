set term post eps color solid enh
set output "try.eps"
set title "Lorenz Equations"
set xlabel "time step"
set ylabel "X"
plot "result.txt" using 1 title "X" with lines


