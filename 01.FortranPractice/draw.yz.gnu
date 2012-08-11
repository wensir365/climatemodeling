set term post eps color solid enh
set output "try.eps"
set title "Lorenz Equations"
set xlabel "Y"
set ylabel "Z"
plot "result.txt" using 2:3 title "Y-Z Phase Space" with lines


