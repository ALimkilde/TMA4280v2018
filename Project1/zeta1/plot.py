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
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Error',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Errors.eps')
plt.figure(2)
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Execution time',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Timings.eps')
plt.show()
