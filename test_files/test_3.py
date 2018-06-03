import numpy as np
import matplotlib.pyplot as plt

# Note higher "kappas" (second arg) result in a _narrower_ distribution

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def dist(x, y):
    """dist"""
    return np.linalg.norm(x-y)

def get_kappa(x, y, nrem, N):
    b = 1.
    d = dist(x,y)
    print "dist = ", d
    lrem = nrem*b
    parameter = np.abs(d-lrem)/b
    if parameter < 2:
        kappa = 30. # really big, have to direct them to each other
    elif parameter < 6:
        kappa = 20
    elif parameter < 10:
        kappa = 10
    elif parameter < 20:
        kappa = 2.
    else:
        kappa = 1.
    return kappa


N = 20
b = 1.
x1 = 0
y1 = 0

x2 = 15
y2 = 0



Thetas = np.zeros(N)
X = np.zeros(N)
Y = np.zeros(N)

for i in range(N/2):
    X[i] = x1
    Y[i] = y1
    X[-(i+1)] = x2
    Y[-(i+1)] = y2

    _, alpha = cart2pol(x2-x1, y2-y1)
    nrem = N - 2*(i+1)
    kappa = get_kappa(np.array([x1,y1]), np.array([x2,y2]),nrem, N)

    alphai_1 = np.random.vonmises(alpha, kappa)
    alphai_2 = np.random.vonmises(np.pi-alpha, kappa)
    Thetas[i] = alphai_1
    Thetas[-(i+1)] = alphai_2


    x1_i = x1 + b * np.cos(alphai_1)
    y1_i = y1 + b * np.sin(alphai_1)

    x2_i = x2 + b * np.cos(alphai_2)
    y2_i = y2 + b * np.sin(alphai_2)

    x1, y1 = x1_i, y1_i
    x2, y2 = x2_i, y2_i

# thetas = np.random.vonmises(np.radians(50), 1.5, 100)
# x, y = np.zeros_like(thetas), np.zeros_like(thetas)
# dx, dy = np.cos(thetas), np.sin(thetas)

# # Convert to x, y, dx, dy...

# for i in range(0,len(thetas)):
#     x[i] = x[i-1] + dx[i]
#     y[i] = y[i-1] + dy[i]


fig, ax = plt.subplots()
ax.plot(X, Y, 'go-')

ax.plot(X[0], Y[0], 'ro')
ax.plot(X[-1], Y[-1], 'ro')
# ax.scatter(x2, y2, 'ro')
ax.quiver(X, Y, np.cos(Thetas), np.sin(Thetas), angles='xy', scale_units='xy', scale=1)
# ax.set(xlim=[-1, 1], ylim=[-1, 1], aspect=1)
ax.axis('off')

plt.show()
