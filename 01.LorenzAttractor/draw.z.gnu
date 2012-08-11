set term post eps color solid enh
set output "try.eps"
set title "Lorenz Equations"
set xlabel "time step"
set ylabel "Z"
plot "result.txt" using 3 title "Z" with lines


