import math
import numpy as np
import random
# def dist_w_pbc(ati, atj, box=[1., 1., 1.]):
#     """
#     calculates the distance between atoms
#     using the pbc
#     input:
#     ati, atj : Atom class instances
#     box: box lengths
#     """
#     dx = atj.x - ati.x
#     dy = atj.y - ati.y
#     dz = atj.z - ati.z
#     dx -= round(dx / box[0]) * box[0]
#     dy -= round(dy / box[1]) * box[1]
#     dz -= round(dz / box[2]) * box[2]
#     return math.sqrt(dx * dx + dy * dy + dz * dz)
#     # min_dists = np.min(np.dstack(((a - b) % bounds, (b - a) % bounds)), axis = 2) dists = np.sqrt(np.sum(min_dists ** 2, axis = 1))

# def dist_wo_pbc(ati, atj, box=[1., 1., 1.]):
#     """
#     calculates the distance between atoms
#     with no pbc
#     input:
#     ati, atj : Atom class instances
#     box: box lengths
#     """
#     dx = atj.x - ati.x
#     dy = atj.y - ati.y
#     dz = atj.z - ati.z
#     return math.sqrt(dx * dx + dy * dy + dz * dz)

def dist_w_pbc(x_1, x_2, bounds=1.1*np.ones(3)):
    """
    calculates the distance between atoms
    with pbc

    input:
    x_1, x_2: np.arrays of 1:3 coordinates
    bounds: np.array 1:3 of the box dimensions
    """
    min_dists = np.min(np.dstack(((x_1 - x_2) % bounds, (x_2 - x_1) % bounds)), axis=2)
    return np.sqrt(np.sum(min_dists ** 2, axis=1))

def dist_wo_pbc(x_1, x_2, bounds=None):
    """
    calculates the distance between atoms
    with no pbc

    input:
    x_1, x_2: np.arrays of 1:3 coordinates
    """
    return np.sqrt(np.sum((x_1-x_2)**2))

def get_pbc(x_1, x_2, bounds):
    """returns the pbc image of the x
    """
    # assert len(x_2) == len(bounds)
    for i in range(3):
        dx = x_1[i] - x_2[i]
        if (dx >  bounds[i] * 0.5):
            x_2[i] += bounds[i]
        if (dx <= -bounds[i] * 0.5):
            x_2[i] -= bounds[i]
    return x_2
def mons_alng_line(x_1, x_2, n, bounds):
    """
    create (n-1) monomers along the line created by ati and atj
    input:
    x_1, x_2 : positions of atoms 1:3 numpy arrays
    n: the length of the polymer chain
    """

    # d = dist_wo_pbc(x_1, x_2, bounds)
    # d_pbc = dist_w_pbc(x_1, x_2, bounds)
    # if d == d_pbc:
    # print "doing mons"
    # print bounds
    x_2 = get_pbc(x_1, x_2, bounds)
    x, dx = np.linspace(x_1[0], x_2[0], n, endpoint=False, retstep=True)
    y, dy = np.linspace(x_1[1], x_2[1], n, endpoint=False, retstep=True)
    z, dz = np.linspace(x_1[2], x_2[2], n, endpoint=False, retstep=True)
    return np.column_stack((x[1:], y[1:], z[1:])), np.linalg.norm([dx, dy, dz])

def three_vector(bond_length, phi, costheta):
    """
    vector from angles
    """
    theta = np.arccos( costheta )
    x = bond_length * np.sin( theta) * np.cos( phi )
    y = bond_length * np.sin( theta) * np.sin( phi )
    z = bond_length * np.cos( theta )
    return (x, y, z)

def random_angles(x1, x2, straight=False):
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    if not straight:
        np.random.seed()
        phi = np.random.uniform(0, np.pi * 2)
        costheta = np.random.uniform(-1, 1)
    if straight:
        r_vec = x1-x2
        r=np.linalg.norm(r_vec)
        costheta = r_vec[2]/r
        phi = np.arctan2(r_vec[1], r_vec[0])
    return phi, costheta




def chain_alng_line(x_1, x_2, n, bounds):
    """
    generates a real chain of a given chain length
    # given two points A and B
    # create random points that are adjacent to A and B
    # interpolate between them repeat until the distance is 1

    """
    b = .8
    N = n
    # flag_verbose = True
    flag_verbose = False

    # x_2 = get_pbc(x_1, x_2, bounds)
    X = np.zeros((N, 3)) # result
    X[0,:] = x_1
    # X[-1,:] = x_2
    x1_i = np.empty_like(x_1)  # our tmp arrays for storing monomers
    x2_i = np.empty_like(x_1)
    dx_1 = np.empty_like(x_1)
    dx_2 = np.empty_like(x_1)

    nleft = N-1
    i = 1
    # print "nleft ", nleft, X
    while nleft >= 1:
        if nleft == 1:
            phi, costheta = random_angles(x_2, x_1, straight=True)
            dx_1 = three_vector(b, phi, costheta)
            X[i, :] = x1_i + dx_1
            nleft -= 1
            break
        # if nleft == 2:
        #     phi, costheta = random_angles(x_1, x_2, straight=True)
        #     dx_1 = three_vector(b, phi, costheta)
        #     phi, costheta = random_angles(x_2, x_1, straight=True)
        #     dx_2 = three_vector(b, phi, costheta)
        #     X[i, :] = x1_i + dx_1
        #     X[-(i+1), :] = x2_i + dx_1
        #     nleft -= 2
        #     break
        # d = dist_w_pbc(x1_i, x2_i, bounds)
        # Nleft = 2 * i + 1
        # dr_left = nleft * b
        phi, costheta = random_angles(x_1, x_2, straight=False)
        dx_1 = three_vector(b, phi, costheta)
        phi, costheta = random_angles(x_1, x_2, straight=False)
        dx_2 = three_vector(b, phi, costheta)
        # if d < dr_left:
        #     # not straight, random
        #     phi, costheta = random_angles(x_1, x_2, straight=False)
        #     dx_1 = three_vector(b, phi, costheta)
        #     phi, costheta = random_angles(x_1, x_2, straight=False)
        #     dx_2 = three_vector(b, phi, costheta)
        # elif d >= dr_left:
        #     phi, costheta = random_angles(x_1, x_2, straight=True)
        #     dx_1 = three_vector(b, phi, costheta)
        #     phi, costheta = random_angles(x_2, x_1, straight=True)
        #     dx_2 = three_vector(b, phi, costheta)
        x1_i = x_1 + dx_1
        x2_i = x_2 + dx_2
        X[i, :] = x1_i
        X[-(i+1), :] = x2_i
        x_1 = x1_i
        x_2 = x2_i
        nleft -= 2
        # print "nleft ", nleft, X
        i += 1
    # return np.column_stack((X[1:], Y[1:], Z[1:])), b
    # print "nleft ", nleft, X
    return X[1:], b
    # return X[1:-2], b

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # for xinit in [0, 2, 3, 5]:
    for xinit in [1]:
        for xend in [9]:
            x_1 = np.array([xinit, 1, 1])
            x_2 = np.array([xend, 1, 1])
            bounds = np.array([10,10,10])
            N = 100
            Coords,_ = chain_alng_line(x_1, x_2, N, bounds)
            # print Coords
            print "chain given N ", N, " result size = ", len(Coords)

            # Coords_mn, _ = mons_alng_line(x_1, x_2, N, bounds)
            # print "mons given N ", N, " result size = ", len(Coords_mn)

            X = Coords[:, 0]
            Y = Coords[:, 1]
            Z = Coords[:, 2]


            # plt.xlim(-10,10)
            # plt.ylim(-10,10)
            # plt.zlim(-10,10)
            ax.set_xlim(-bounds[0], bounds[0])
            ax.set_ylim(-bounds[0], bounds[0])
            ax.set_zlim(-bounds[0], bounds[0])
            ax.scatter(X[0], Y[0], Z[0], 'ro',color='red',marker='D',s=50)
            ax.scatter(X[-1], Y[-1], Z[-1], 'ro',color='red',marker='D',s=50)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.plot(X, Y, Z, 'bo-')
    plt.show()

