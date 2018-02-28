# Molecular Dynamics Simulation Program [Project 2018 EIA (II)]

## Getting Started

### Prerequisites
This program is written in Fortran 90 and the visualizations tools in Python 2.4.  You need to have a fortran compiler, **gcc** works fine and a fesh installation of Python 2.4.

The following Python libraries are required to plot the data:

```
matplotlib numpy
```

### How to compile it?
A Makefile has been created to make the work easier. When the make is called:
* The modules are going to be compiled as fortran objects: 
```
$gfortran -c * .f90
```
* The program will be compiled as: 
```
$gfortran -c md * .o -O3
```
which means that an optimization of level 3 will be performed. *
* To execute the program, the file with the parameters is going to be called in the command. ./md param.dat
* Once the results are stored in different files, plots will be created using the python and gnuplot scripts.
* Finally, all the results will be printed on screen.

## Running the program


