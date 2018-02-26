echo Generando script gnuplot para $1
sed "s/TEMPLATE/$1/" TEMPLATE.gnuplot > $1.gnuplot
