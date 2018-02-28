# Molecular Dynamics Simulation Program [Project 2018 EIA (II)]

## Getting Started

### Prerequisites
This program is written in Fortran 90 and the visualizations tools in Python 2.4.  You need to have a fortran compiler, (**gcc** should be fine) and a fresh installation of Python 2.4.

The following Python libraries are required to plot the data:

```
matplotlib numpy
```

Also you need **make** to execute correctly the makefile. **make** package is included in most Linux distributions, if you need to download it run the following commands depending on your system:

**Fedora 27** 
```
sudo dnf install make
```
**Debian Stretch/Ubuntu 16.04**
```
sudo apt-get install make
```
**Arch/Manjaro**
```
sudo pacman -S make
```
### How to download the code
If you have an installation of git in your computer you can clone this repository using the clone command:
```
git clone https://github.com/EIA-Master/Project_2018_II
```
We are working now with the **Tidy** branch, so is needed to change the branch to work with the proper files:
You can also download a .zip with the code of this Branch using the GitHub download button.
### How to compile it?
A Makefile has been created to make the work easier. When the make is called using the command:
```
make
```
* The modules are going to be compiled as fortran objects: 
```
$gfortran -c * .f90
```
* The program will be compiled as: 
```
$gfortran -c md * .o -O3
```
which means that an optimization of level 3 will be performed.
* To execute the program, the file with the parameters is going to be called in the command. ./md param.dat
* Once the results are stored in different files, plots will be created using the python and gnuplot scripts.
* Finally, all the results will be printed on screen.

## Running the program

### Parameters

There is a file called param.dat, this file contains all the variables that you can change to perform the dynamic. You must introduce the variables in reduced units.
Only Verlet (**Ver**) and Euler (**Eul**) integrators are implemented.

### Execution

It's recommened to execute the file with the command:
```
make datum
```
But you also can run the code running the **md** binary.

### Plotting
If you want to plot the output you need to execute the following command:
```
make plot
```
To represent the dynamic you need a visualization program. For example:
* **VMD:** http://www.ks.uiuc.edu/Research/vmd/
* **wxMacMolPlt**: https://brettbode.github.io/wxmacmolplt/downloads.html
The xyz coordinates are stored at **traj.xyz** file.

## Authors
* Alejandro C. - *Initial geometry and makefile*
* Cristian P. - *Integrators and Andersen implementation*
* Erika F. - *GDR implementation and PBC*
* Zazo M. - *Statistical analysis and plotting*
* Sergio G. - *GitHub management, Main program and Lennard Jones implementation.*
