import matplotlib.pyplot as plt
import numpy as np



time = []
Ptot = []
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
            Epot.append(words[2])
            Ekin.append(words[3])
            Etot.append(words[4])
            
        except:
            pass

x_list = time

plt.plot(x_list, Epot)
plt.title('test visualisation of data')
plt.xlabel('time (unit?)')
plt.ylabel('Energy (au)')
plt.show()
