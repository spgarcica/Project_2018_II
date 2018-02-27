# Zazo Meijs

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def Round2(value):
    return str(round(value,2))[0:4]


velocities = []
n = 0
with open('testvel', 'r') as file:

    for line in file:
        n += 1


        values = line.split()
        velocities.append([float(values[i]) for i in range(len(values))])
                          

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

x_velocities = [velocities[i][0] for i in range(len(velocities))]

# the histogram of the data
n, bins, patches = plt.hist(x_velocities, 50,normed=1, facecolor='green', alpha=0.75, label='Results')
(mu, sigma) = norm.fit(x_velocities)

# add a 'best fit' line
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=1, label = 'fit mu:%s, sigma:%s'%(Round2(mu), Round2(sigma)))

plt.xlabel('Velocities(au)')
plt.ylabel('Percentage()')
plt.title('Normalized distribution of velocities in the x direction')
plt.legend()
plt.grid(True)
plt.savefig('velocity_distribution.png', bbox_inches='tight')

plt.show()
