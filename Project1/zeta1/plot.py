import matplotlib.pyplot as plt
import numpy as np


Nvec = np.power(2,[1,2,3,4,5,6])

for n in Nvec:
    E = np.loadtxt('res/Errors%0.0f.txt' %n)
    T = np.loadtxt('res/Times%0.0f.txt' %n)
    plt.figure(1)
    plt.loglog(E[:,0],E[:,1],lw = 2,label = 'Np = %0.0f' %n)

    plt.figure(2)
    plt.loglog(T[:,0],T[:,1],lw = 2,label = 'Np = %0.0f' %n)

plt.figure(1)
plt.legend(loc = 'best')
plt.ylabel('Error')
plt.xlabel('n')
plt.figure(2)
plt.legend(loc = 'best')
plt.ylabel('Execution time')
plt.xlabel('n')
plt.show()
