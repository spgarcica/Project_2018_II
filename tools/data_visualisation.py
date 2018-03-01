import matplotlib.pyplot as plt
import numpy as np




def Bin_Analysis( array):
    '''
    The analasys takes the last 40% of the data where from visual observation
    we see that it has equibrilated. And split in 40 blocks of which 30 are used
    for binning statistics and then are left unused to decorelate the bins.
    @In : The complete array
    @out: Average, Stdv (Standart deviation)
    '''

    if type(array) != np.ndarray:   # If it's not a numpy array
        array = np.array(array)     # Make it a Numpy array ;)

    averages = []
    N = len(array)  
    for n in range(6*N/10, (N+1)-N/25, N/25): #Remember that in Python 2.7 int/int=int 
        average_n = np.average(array[n:(n+3*N/100)])
        averages.append(average_n)

    Average = np.average(averages)
    Stdv = 2*np.std(averages)/np.sqrt(len(averages)) # 95% certainty interval

    #print Average, Stdv
    return int(np.round(Average)), int(np.round(Stdv))

                            
    

# Read in Energies and Momenta

time = []   # time in units of frame
Ptot = []   # total Momentum
Press= []   # Pressure
Temp = []   # Temperature
Epot = []   # potential Energy
Ekin = []   # kinetic Energy
Etot = []   # total Energy

with open("ener.xyz", 'r') as data:
    for line in data:
        words = line.split()
        try:
            words = [float(words[i]) for i in range(len(words))]
            time.append(words[0])
            Ptot.append(words[1])
            Press.append(words[2])
            Epot.append(words[3])
            Ekin.append(words[4])
            Etot.append(words[5])
            
        except:
            pass

with open("temperature.data", 'r') as data:
    for line in data:
        words = line.split()
        try:
            words = [float(words[i]) for i in range(len(words))]
            Temp.append(words[0])
            
        except:
            pass


# Plot Energies
x_list = time
plt.plot(x_list, Etot, label="Total Energy: %i (%i))"%(Bin_Analysis(Etot)))
plt.plot(x_list, Ekin, label="Kinectic Energy: %i (%i))"%(Bin_Analysis(Ekin)))
plt.plot(x_list, Epot, label="Potential Energy: %i (%i))"%(Bin_Analysis(Epot)))
plt.legend()

plt.title('Visualisation of Energies, values given are average of last 40% with 95% certainty interval')
plt.xlabel('time (freme)')
plt.ylabel('Energy (reduced units)')

plt.savefig('energies.png', bbox_inches='tight')
#plt.show()
plt.clf()

#'''
# Plot Temperature
x_list = time
plt.plot(x_list, Temp)

plt.title('Temperature, equalizes at: %i (%i)'%(Bin_Analysis(Temp)))
plt.xlabel('time (frame)')
plt.ylabel('Temperature (R.U)')

plt.savefig('temperature.png', bbox_inches='tight')
#plt.show()
plt.clf()
#'''

#'''
# Plot Pressure
x_list = time
plt.plot(x_list, Press)

plt.title('Pressure, equalizes at: %i (%i)'%(Bin_Analysis(Press)))
plt.xlabel('time (frame)')
plt.ylabel('Pressure (au)')

plt.savefig('pressure.png', bbox_inches='tight')
#plt.show()
plt.clf()
#'''
