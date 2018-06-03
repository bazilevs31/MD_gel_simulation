import random
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

import numpy as np
import math

def dist(x1, y1, x2, y2):
    """Summary
    calculates distance between two points
    Args:
    x1 (float): point
    y1 (float): point
    x2 (float): point
    y2 (float): point

    Returns:
    float: distance
    """
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)


N = 50
b = 1.
x1 = 0
y1 = 0

x2 = 15
y2 = 0

plt.xlim(-30,30)
plt.ylim(-30,30)

if N*b <= dist(x1, y1, x2, y2):
    print "error should be longer"
else:
    X = np.empty(N)
    Y = np.empty(N)


plt.plot(x1, y1, 'ro',markersize=20)
plt.plot(x2, y2, 'ro',markersize=20)
for i in range(N/2):
    X[i] = x1
    Y[i] = y1
    X[-(i+1)] = x2
    Y[-(i+1)] = y2

    alphai_1 = math.pi*random.random()
    alphai_2 = math.pi*random.random()

    x1_i = x1 + b * math.cos(alphai_1)
    y1_i = y1 + b * math.sin(alphai_1)

    x2_i = x2 + b * math.cos(alphai_2)
    y2_i = y2 + b * math.sin(alphai_2)

    d = dist(x1_i, y1_i, x2_i, y2_i)
    Nleft = N - 2 * (i+1)
    dr_left = Nleft * b
    if d < dr_left:
        print "normal step, dist", d, " dist left = ", dr_left, " 1:x,y ", x1_i, y1_i, " 2:x,y ", x2_i, y2_i
        x1, y1 = x1_i, y1_i
        x2, y2 = x2_i, y2_i
    elif np.allclose(d, dr_left):
        # connect with a straight line
        print "straight line step, dist", d, " dist left = ", dr_left
        x1, y1 = x1_i, y1_i
        x2, y2 = x2_i, y2_i
        for j in range(Nleft):
            x1_i = x1 + b * i
            y1_i = y1 + b * i
            X[i+j+1], Y[i+j+1] = x1_i, y1_i
            print "steps ", x1_i, y1_i
        break
    elif d > dr_left:
        # do some magic with reiterating the last step
        print "too close line step, dist", d, " dist left = ", dr_left
        # x1, y1 = x1_i, y1_i
        # x2, y2 = x2_i, y2_i
        x = np.linspace(x1, x2, Nleft, endpoint=False, retstep=False)
        y = np.linspace(y1, y2, Nleft, endpoint=False, retstep=False)
        for j in range(Nleft):
            X[i+j+1], Y[i+j+1] = x[j], y[j]
            print "steps ", x[j], y[j]
        break
    else:
        print "error"


print X
print Y
plt.plot(X, Y, 'bo-')

plt.show()
