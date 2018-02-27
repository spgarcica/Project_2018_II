# Project_2018_II

The code presented in this repository performs a molecular dynamics simulation.
In this simulation, two integrator algorithms can be choose, Euler or Velocity Verlet and the temperature is controlled using the Andersen thermostat. There is no pressure coupling. 

## How does it work?
[Esquema de como lo hace]

## The code
The main program needs the following modules:
* **mtmod**: Includes a random number generator. The generation of random numbers is needed for the Andersen thermostat and for generating random initial velocities.
* **normaldist**:
* **pbcmod**: Contains the subroutines that are necessary to implement the periodic boundary conditions of the system.
* **rdfmod**: Calculates the radial distribution function.
* **LJ_forcemod**: Generates the matrix containing the forces that act in every atom and the potential energy.
* **eulermod**: Includes the Euler integrator algorithm.
* **velocity_verletmod**: Include the Velocity-Verlet integrator algorithm.
* **andersenmod**: Includes the thermostat.
* **gaussmod**: In order to check if the velocities obtained follow a normal distribution, this module generates a gaussian curve to make the comparisson with the obatained data.
* **geommod**: Generates a cubic simple crystalline structure as an initial geometry. This structure will be melted as to have a fluid.
* **initial_velocitymod**: Generates random but normally distributed velocities, around the value of sigma.
* **momentumod**:
* **repairmod**:
* **statisticalsmod**:
* **stmod**:
