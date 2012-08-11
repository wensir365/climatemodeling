set term post eps color solid enh
set output "try.eps"
set title "Lorenz Equations"
splot "result.txt" title "X-Y-Z Phase Space" with lines


