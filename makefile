CODE = $(shell find -name '*.f90')
MTMOD = $(shell find -name '*mtmod.f90')
NORMALDIST = $(shell find -name '*normaldist.f90')
MOMENTUMOD = $(shell find -name '*momentumod.f90')
I_VEL = $(shell find -name '*initial_velocitymod.f90')
ANDERSEN = $(shell find -name '*andersenmod.f90')
REPAIR = $(shell find -name '*repairmod.f90')
GEOMMOD = $(shell find -name '*geommod.f90')
PBCMOD = $(shell find -name '*pbcmod.f90')
GAUSS = $(shell find -name '*gaussmod.f90')
STD = $(shell find -name '*stdmod.f90')
STATISTIC = $(shell find -name '*statisticalsmod.f90')
RDF = $(shell find -name '*rdfmod.f90')
LJ = $(shell find -name '*LJ_forcemod.f90')
VV = $(shell find -name '*velocity_verlet.f90')
EULER = $(shell find -name '*eulermod.f90')
MAIN = $(shell find -name '*md.f90')
PARAM = $(shell find -name '*param.dat')
DATA = ener.xyz gauss.dat g_function.txt temperature.data testvel traj.xyz
DATA_SCRIPT_PYTHON = ener.xyz
DATA_SCRIPT_PYTHON_DIST_VEL = testvel
DATA_SCRIPT_GNUPLOT = g_function.txt
SCRIPT_PYTHON = $(shell find -name '*data_visualisation.py')
SCRIPT_PYTHON_DIST_VEL = $(shell find -name '*velocities_distribution.py')
SCRIPT_GNUPLOT = $(shell find -name '*.gnuplot')
FIGURE_SCRIPT_PYTHON = energies.png temperature.png pressure.png
FIGURE_SCRIPT_PYTHON_DIST_VEL = velocity_distribution.png
FIGURE_SCRIPT_GNUPLOT = RDF.eps
TARGET = md


$(TARGET) : $(TARGET).o
	gfortran -o $(TARGET) *.o -O3

$(TARGET).o : $(CODE)
	gfortran -c $(MTMOD)
	gfortran -c $(NORMALDIST)
	gfortran -c $(MOMENTUMOD)
	gfortran -c $(I_VEL)
	gfortran -c $(ANDERSEN)
	gfortran -c $(REPAIR)
	gfortran -c $(GEOMMOD)
	gfortran -c $(PBCMOD)
	gfortran -c $(GAUSS)
	gfortran -c $(STD)
	gfortran -c $(STATISTIC)
	gfortran -c $(RDF)
	gfortran -c $(LJ)
	gfortran -c $(VV)
	gfortran -c $(EULER)
	gfortran -c $(MAIN)

$(DATA) : $(TARGET) $(PARAM)
	./$(TARGET)

$(FIGURE_SCRIPT_PYTHON) : $(DATA_SCRIPT_PYTHON) $(SCRIPT_PYTHON)
	python2 $(SCRIPT_PYTHON)

$(FIGURE_SCRIPT_PYTHON_DIST_VEL) : $(DATA_SCRIPT_PYTHON_DIST_VEL) $(SCRIPT_PYTHON_DIST_VEL)
	python2 $(SCRIPT_PYTHON_DIST_VEL)

$(FIGURE_SCRIPT_GNUPLOT) : $(DATA_SCRIPT_GNUPLOT) $(SCRIPT_GNUPLOT)
	gnuplot $(SCRIPT_GNUPLOT)

## datum : generate data files about MD simulation
.PHONY : datum
datum : $(DATA)

## plot : generate figures about various magnitudes of interest in MD
.PHONY : plot
plot : $(FIGURE_SCRIPT_PYTHON) $(FIGURE_SCRIPT_PYTHON_DIST_VEL) $(FIGURE_SCRIPT_GNUPLOT)

## compilation : compilation of the program
.PHONY : compilation
compilation : $(TARGET)

## variables : print variables
.PHONY : variables
variables :
	@echo CODE:$(CODE)
	@echo MTMOD:$(MTMOD)
	@echo NORMALDIST:$(NORMALDIST)
	@echo MOMENTUMOD:$(MOMENTUMOD)
	@echo I_VEL:$(I_VEL)
	@echo ANDERSEN:$(ANDERSEN)
	@echo REPAIR:$(REPAIR)
	@echo GEOMMOD:$(GEOMMOD)
	@echo PBCMOD:$(PBCMOD)
	@echo GAUSS:$(GAUSS)
	@echo STD:$(STD)
	@echo STATISTIC:$(STATISTIC)
	@echo RDF:$(RDF)
	@echo LJ:$(LJ)
	@echo VV:$(VV)
	@echo EULER:$(EULER)
	@echo MAIN:$(MAIN)
	@echo PARAM:$(PARAM)
	@echo DATA:$(DATA)
	@echo DATA_SCRIPT_PYTHON:$(DATA_SCRIPT_PYTHON)
	@echo DATA_SCRIPT_PYTHON_DIST_VEL:$(DATA_SCRIPT_PYTHON_DIST_VEL)
	@echo DATA_SCRIPT_GNUPLOT:$(DATA_SCRIPT_GNUPLOT)
	@echo SCRIPT_PYTHON:$(SCRIPT_PYTHON)
	@echo SCRIPT_PYTHON_DIST_VEL:$(SCRIPT_PYTHON_DIST_VEL)
	@echo SCRIPT_GNUPLOT:$(SCRIPT_GNUPLOT)
	@echo FIGURE_SCRIPT_PYTHON:$(FIGURE_SCRIPT_PYTHON)
	@echo FIGURE_SCRIPT_PYTHON_DIST_VEL:$(FIGURE_SCRIPT_PYTHON_DIST_VEL)
	@echo FIGURE_SCRIPT_GNUPLOT:$(FIGURE_SCRIPT_GNUPLOT)
	@echo TARGET:$(TARGET)

## clean : remove auto-generated files
.PHONY : clean
clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(DATA)
	rm -f $(FIGURE_SCRIPT_PYTHON)
	rm -f $(FIGURE_SCRIPT_PYTHON_DIST_VEL)
	rm -f $(FIGURE_SCRIPT_GNUPLOT)
	rm -f $(TARGET)

## help : provide some instructions useful to use the makefile
.PHONY : help
help :
	@sed -n 's/^##//p' makefile
