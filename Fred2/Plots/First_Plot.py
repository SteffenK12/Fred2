from copy import deepcopy

__author__ = 'Steffen'
import matplotlib.pyplot as plt
import math
import numpy as np

def f(x,y):
    return x

x = np.array([1, 2, 3, 4, 5, 6, 7])
y = np.array([2, 4, 7, 7, 7, 4, 2])
#X, Y = np.meshgrid(x,y)
#Z = f(x,y)




low = min(y)
high = max(y)

colormap = plt.get_cmap('Reds')
fig = plt.figure()
ax1 = plt.subplot2grid((1, 1), (0, 0))

#plt.pcolormesh(X, Y, Z, cmap=colormap)
ax1.plot(x, y, "-", label='random data')
ax1.grid(True)
#plt.xlabel('x')
#plt.ylabel('y')
plt.ylim([0, math.ceil(high + 0.5 * (high - low))])
plt.subplots_adjust(top=0.95)
plt.title('Test Graph')
#plt.legend()
plt.show()

