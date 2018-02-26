set title "Radial Distribution Function"
set xlabel "Distance"
set ylabel "g(r)"
unset key

set terminal postscript eps color; set o 'RDF.eps'

plot 'g_function.txt' w l 
