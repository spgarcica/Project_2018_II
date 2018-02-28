# Molecular Dynamics Simulation Program [Project 2018 EIA (II)]

## How to compile it?
A Makefile has been created to make the work easier. When the make is called:
* The modules are going to be compiled as fortran objects: $gfortran -c * .f90
* The program will be compiled as: $gfortran -c md * .o -O3, which means that an optimization of level 3 will be performed.
* To execute the program, the file with the parameters is going to be called in the command. ./md param.dat
* Once the results are stored in different files, plots will be created using the python and gnuplot scripts.
* Finally, all the results will be printed on screen.
