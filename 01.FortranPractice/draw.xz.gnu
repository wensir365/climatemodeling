set term post eps color solid enh
set output "try.eps"
set title "Lorenz Equations"
set xlabel "X"
set ylabel "Z"
plot "result.txt" using 1:3 title "X-Z Phase Space" with lines


