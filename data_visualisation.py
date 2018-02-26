import matplotlib.pyplot as plt
import numpy as np



time = []
Ptot = []
Press= []
Epot = []
Ekin = []
Etot = []

with open("temp_res.xyz", 'r') as data:
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

x_list = time

plt.plot(x_list, Etot, label="Total Energy")
plt.plot(x_list, Ekin, label="Kinectic Energy")
plt.plot(x_list, Epot, label="Potential Energy")
plt.legend()

plt.title('Test visualisation of Energies')
plt.xlabel('time (unit?)')
plt.ylabel('Energy (au)')
plt.show()
