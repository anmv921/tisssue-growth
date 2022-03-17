set terminal qt
reset
plot "props.txt" u 1:3 w linespoints pt 6
set xlabel "t"
set ylabel "P"
unset key
replot