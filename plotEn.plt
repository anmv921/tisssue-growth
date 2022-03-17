set terminal qt
reset
plot "props.txt" u 1:2 w linespoints pt 6
set xlabel "t"
set ylabel "E"
unset key
replot