set xlabel 'TEMPLATE'
set ylabel 'Frame'
set terminal png
set output 'TEMPLATE.png'
plot 'TEMPLATE.dat' u 1:2 w l
