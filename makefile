CODE = $(shell find -name '*.f90')
OBJECTS = $(CODE:.f90=.o)
DATA = Ek.dat
DATA_SUFFIX = $(DATA:.dat=) 
GNUPLOT_DATA = $(DATA:.dat=.gnuplot)
TARGET = MainProgMod


$(TARGET) : $(TARGET).o
	gfortran -o $(TARGET) $(CODE)

$(TARGET).o : $(CODE)
	gfortran -c $(CODE)

$(DATA) : 
	./$(TARGET)

$(GNUPLOT_DATA) : $(DATA)
	./genera_gnuplot.sh $(DATA_SUFFIX)
	gnuplot $@

## datum : generate data files about MD simulation
.PHONY : datum
datum : $(DATA)

## plot : generate files about various magnitudes of interest in MD
.PHONY : plot
plot : $(GNUPLOT_DATA)

## compilation : compilation of the program
.PHONY : compilation
compilation : $(TARGET)

## variables : print variables
.PHONY : variables
variables :
	@echo CODE:$(CODE)
	@echo OBJECTS:$(OBJECTS)
	@echo DATA:$(DATA)
	@echo DATA_SUFFIX:$(DATA_SUFFIX)
	@echo TARGET:$(TARGET)

.PHONY : help
help :
	@sed -n 's/^##//p' makefile
