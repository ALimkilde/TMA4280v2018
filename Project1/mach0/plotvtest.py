import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt('vtest.txt')

plt.figure()
plt.loglog(M[:,0],M[:,1])
plt.show()
