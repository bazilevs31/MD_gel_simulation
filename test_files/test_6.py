import numpy as np
import matplotlib.pyplot as plt

# Note higher "kappas" (second arg) result in a _narrower_ distribution

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def cart2spher(x, y, z):
    rho = np.sqrt(x**2 + y**2 + z**2)
    theta = z/float(rho)
    phi = np.arctan2(y, x)
    return(rho, phi, theta)


def dist_3d(x1, y1,z1, x2, y2,z2):
    """Summary

    Args:
        x1 (float): point
        y1 (float): point
        z1 (float): point
        x2 (float): point
        y2 (float): point
        z2 (float): point

    Returns:
        float: dist
    """
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)




def dist(x, y):
    """dist"""
    return np.linalg.norm(x-y)

def get_kappa(x, y, nrem, N):
    b = 1.
    d = dist(x,y)
    # print "dist = ", d
    lrem = nrem*b
    parameter = np.abs(d-lrem)/b
    # if parameter > 10:
    #     kappa = 0.1
    # else:
    #     kappa = 20
    if parameter < 2:
        kappa = 20. # really big, have to direct them to each other
    elif parameter < 6:
        kappa = 16
    elif parameter < 10:
        kappa = 5
    elif parameter < 15:
        kappa = .9
    else:
        kappa = 0.1
    return kappa


N = 100
b = 1.
x1 = 0
y1 = 0
z1 = 0

x2 = 15
y2 = 0
z2 = 0


Thetas = np.zeros(N)
X = np.zeros(N)
Y = np.zeros(N)
Z = np.zeros(N)
X[0] = x1
Y[0] = y1
Z[0] = z1

X[-1] = x2
Y[-1] = y2
Z[-1] = z2

# for N in [20, 25, 30, 40, 50]:
# for N in [50]:
for nsteps in range(1):
    for i in range(1, N-1):
        _, phi_direction, theta_direction = cart2spher(x2-x1, y2-y1, z2-z1)
        nrem = N - (i+2)
        kappa = get_kappa(np.array([x1,y1]), np.array([x2,y2]),nrem, N)

        alphai_1 = np.random.vonmises(phi_direction, kappa)
        thetai_1 = np.random.vonmises(theta_direction, kappa)
        Thetas[i] = alphai_1


        x1_i = x1 + b * np.cos(alphai_1)
        y1_i = y1 + b * np.sin(alphai_1)

        x1_i = x1 + b * np.cos(alphai_1) * np.cos(thetai_1)
        y1_i = y1 + b * np.sin(alphai_1) * np.cos(thetai_1)
        z1_i = z1 + b * np.sin(thetai_1)

        print "nrem = ", nrem, "dist = ", dist_3d(x1_i, y1_i, z1_i, x1, y1, z1)


        x1, y1,z1 = x1_i, y1_i,z1_i
        X[i] = x1
        Y[i] = y1
        Z[i] = z1
    print "nrem = ", nrem, "dist = ", dist_3d(x1_i, y1_i, z1_i, x2, y2, z2)

        # print "N", N"nrem = ", nrem, "dist = ", dist_2d(x1_i, y1_i, x2, y2)
        # print "N ", N, "dist = ", dist_3d(x1_i, y1_i, x2, y2)


# thetas = np.random.vonmises(np.radians(50), 1.5, 100)
# x, y = np.zeros_like(thetas), np.zeros_like(thetas)
# dx, dy = np.cos(thetas), np.sin(thetas)

# # Convert to x, y, dx, dy...

# for i in range(0,len(thetas)):
#     x[i] = x[i-1] + dx[i]
#     y[i] = y[i-1] + dy[i]
def writexyz(X, Y, Z, filename='mol_default.xyz'):
        outfile = filename
        nparticles = len(X)
        with open(outfile, 'w') as f:
            f.write(str(nparticles) + '\n')
            f.write('elementary file' + '\n')
            f.write("O {x} {y} {z}\n".format(x=X[0], y=Y[0], z=Z[0]))
            for i in range(1, nparticles-1):
                # i += 1
                f.write("C {x} {y} {z}\n".format(x=X[i], y=Y[i], z=Z[i]))
            f.write("O {x} {y} {z}\n".format(x=X[-1], y=Y[-1], z=Z[-1]))
        return None


writexyz(X, Y, Z)
# from mpl_toolkits.mplot3d import proj3d
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# # ax.scatter(x1[0],x1[1],x1[2], color='red',marker='D',s=50)
# # ax.scatter(x2[0],x2[1],x2[2], color='red',marker='D',s=50)
# ax.plot(X, Y, Z, 'go-')
# ax.scatter(X[0], Y[0], Z[0], 'rD')
# ax.scatter(X[-1], Y[-1], Z[-1], 'D')

# # ax.plot(X[0], Y[0], 'ro')
# # ax.plot(X[-1], Y[-1], 'ro')
# # # ax.scatter(x2, y2, 'ro')
# # ax.quiver(X, Y, np.cos(Thetas), np.sin(Thetas), angles='xy', scale_units='xy', scale=1)
# # # ax.set(xlim=[-1, 1], ylim=[-1, 1], aspect=1)
# # ax.axis('off')

# plt.show()
# fig, ax = plt.subplots()
# ax.plot(X, Y, 'go-')

# ax.plot(X[0], Y[0], 'ro')
# ax.plot(X[-1], Y[-1], 'ro')
# # # ax.scatter(x2, y2, 'ro')
# # ax.quiver(X, Y, np.cos(Thetas), np.sin(Thetas), angles='xy', scale_units='xy', scale=1)
# # # ax.set(xlim=[-1, 1], ylim=[-1, 1], aspect=1)
# # ax.axis('off')

# plt.show()
