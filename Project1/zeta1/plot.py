import matplotlib.pyplot as plt
import numpy as np


Nvec = np.power(2,[0,1,2,3,4])

for n in Nvec:
    E = np.loadtxt('res/Errors%0.0f.txt' %n)
    T = np.loadtxt('res/Times%0.0f.txt' %n)
    plt.figure(1)
    plt.loglog(E[:,0],E[:,1],lw = 2,label = 'Np = %0.0f' %n)

    plt.figure(2)
    plt.plot(T[:,0],T[:,1],lw = 2,label = 'Np = %0.0f' %n)

plt.figure(1)
plt.loglog(E[0:10,0],0.4/E[0:10,0],'k--',lw = 2,label = 'O(1/n)')
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Error',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Errors.eps')
plt.figure(2)
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Execution time',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Timings.eps')

Tvec = np.zeros(len(Nvec))
for i,n in enumerate(Nvec):
    T = np.loadtxt('res/Times%0.0f.txt' %n)
    Tvec[i] =T[7,1]

plt.figure(3)
plt.plot(Nvec,Tvec[0]/Tvec)
#plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Execution time',fontsize = 12)
plt.xlabel('nprocs',fontsize = 12)
plt.savefig('SpeedUp.eps')


for n in Nvec:
    E = np.loadtxt('Lille/Errors%0.0f.txt' %n)
    T = np.loadtxt('Lille/Times%0.0f.txt' %n)
    plt.figure(3)
    plt.loglog(E[:,0],E[:,1],lw = 2,label = 'Np = %0.0f' %n)

    plt.figure(4)
    plt.plot(T[:,0],T[:,1],lw = 2,label = 'Np = %0.0f' %n)

plt.figure(3)
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Error',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Errors_Lille.eps')
plt.figure(4)
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Execution time',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Timings_Lille.eps')
