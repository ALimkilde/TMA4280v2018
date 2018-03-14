import matplotlib.pyplot as plt
import numpy as np


Nvec = np.power(2,[0,1])
Mvec = np.power(2,[0,1,2])

marks = ['o-','s--','^-.','-.']
colors = ['k','r','g','b','y','m']
i = 0
j = 0

for m in Mvec:
    j = 0
    for n in Nvec:
        E = np.loadtxt('res2/Errors_np%0.0f_m%0.0f.txt' %(n,m))
        T = np.loadtxt('res2/Times_np%0.0f_m%0.0f.txt' %(n,m))
        #plt.figure(1)
        #plt.loglog(E[:,0],E[:,1],marks[i],lw = 2,label = 'Np = %0.0f' %n)

        plt.figure(2)
        plt.plot(T[:,0],T[:,1],marks[j],color = colors[i],lw = 2,label = 'Np = %0.0f, m = %0.0f' %(n,m))
        j += 1
    i += 1
plt.figure(2)
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Execution time',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Timings.eps')

Tvec = np.zeros(len(Mvec))
i=0
for m in Mvec:
    j = 0
    for n in Nvec:
        T = np.loadtxt('res2/Times_np%0.0f_m%0.0f.txt' %(n,m))
        #plt.figure(1)
        #plt.loglog(E[:,0],E[:,1],marks[i],lw = 2,label = 'Np = %0.0f' %n)

        Tvec[i] = T[-1,1]
        j += 1
    i += 1
plt.figure(3)
plt.plot(Mvec,Tvec)
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Execution time',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Timings2.eps')
plt.show()

i = 0;
for m in Mvec:
    j = 0
    for n in Nvec:
        E = np.loadtxt('Lille/Errors_np%0.0f_m%0.0f.txt' %(n,m))
        T = np.loadtxt('Lille/Times_np%0.0f_m%0.0f.txt' %(n,m))
        #plt.figure(1)
        #plt.loglog(E[:,0],E[:,1],marks[i],lw = 2,label = 'Np = %0.0f' %n)

        plt.figure(3)
        plt.plot(T[:,0],T[:,1],marks[j],color = colors[i],lw = 2,label = 'Np = %0.0f, m = %0.0f' %(n,m))
        j += 1
    i += 1
plt.figure(3)
plt.legend(loc = 'best',fontsize = 12)
plt.ylabel('Execution time',fontsize = 12)
plt.xlabel('n',fontsize = 12)
plt.savefig('Timings_Lille.eps')
