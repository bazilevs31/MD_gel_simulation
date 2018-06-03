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
    # print( "doing mons")
    # print( bounds)
    x_2 = get_pbc(x_1, x_2, bounds)
    x, dx = np.linspace(x_1[0], x_2[0], n, endpoint=False, retstep=True)
    y, dy = np.linspace(x_1[1], x_2[1], n, endpoint=False, retstep=True)
    z, dz = np.linspace(x_1[2], x_2[2], n, endpoint=False, retstep=True)
    return np.column_stack((x[1:], y[1:], z[1:])), np.linalg.norm([dx, dy, dz])

def tether_alng_line(x_1, x_2, n, bounds):
    """
    create (n-1) monomers along the line created by ati and atj
    input:
    x_1, x_2 : positions of atoms 1:3 numpy arrays
    n: the length of the polymer chain
    """

    # d = dist_wo_pbc(x_1, x_2, bounds)
    # d_pbc = dist_w_pbc(x_1, x_2, bounds)
    # if d == d_pbc:
    # print( "doing mons")
    # print( bounds)
    # x_2 = get_pbc(x_1, x_2, bounds)
    x, dx = np.linspace(x_1[0], x_2[0], n, endpoint=False, retstep=True)
    y, dy = np.linspace(x_1[1], x_2[1], n, endpoint=False, retstep=True)
    z, dz = np.linspace(x_1[2], x_2[2], n, endpoint=False, retstep=True)
    return np.column_stack((x[1:], y[1:], z[1:])), np.linalg.norm([dx, dy, dz])




def random_three_vector(bond_length):
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    np.random.seed()
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = bond_length * np.sin( theta) * np.cos( phi )
    y = bond_length * np.sin( theta) * np.sin( phi )
    z = bond_length * np.cos( theta )
    return (x, y, z)


def chain_alng_line(x_1, x_2, n, bounds):
    """
    generates a real chain of a given chain length
    # given two points A and B
    # create random points that are adjacent to A and B
    # interpolate between them repeat until the distance is 1

    """
    b = .8
    N = n
    x_2 = get_pbc(x_1, x_2, bounds)
    xx1, yy1, zz1 = x_1[0], x_1[1], x_1[2]
    xx2, yy2, zz2 = x_2[0], x_2[1], x_2[2]
    # flag_verbose = True
    flag_verbose = False

    def dist(x1, y1, z1, x2, y2, z2):
        return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)


    if N*b <= dist(xx1, yy1, zz1, xx2, yy2, zz2):
        # print( "error should be longer")
        raise ValueError("error should be longer")
    else:
        X = np.zeros(N)
        Y = np.zeros(N)
        Z = np.zeros(N)

    X[0] = xx1
    Y[0] = yy1
    Z[0] = zz1
    for i in range(N/2):
        (deltax, deltay, deltaz) = random_three_vector(b)
        (xx1_i, yy1_i, zz1_i) = (xx1 + deltax, yy1 + deltay, zz1 + deltaz)
        (deltax, deltay, deltaz) = random_three_vector(b)
        (xx2_i, yy2_i, zz2_i) = (xx2 + deltax, yy2 + deltay, zz2 + deltaz)

        d = dist(xx1_i, yy1_i, zz1_i, xx2_i, yy2_i, zz2_i)
        Nleft = N - 2 * (i+1)
        dr_left = Nleft * b
        if d < dr_left:
            if flag_verbose:
                print( "normal step, dist", d, " dist left = ", dr_left, " 1:x,y ", xx1_i, yy1_i, " 2:x,y ", xx2_i, yy2_i)

            xx1, yy1, zz1 = xx1_i, yy1_i, zz1_i
            xx2, yy2, zz2 = xx2_i, yy2_i, zz2_i
            X[i + 1] = xx1
            Y[i + 1] = yy1
            Z[i + 1] = zz1
            X[-(i + 1)] = xx2
            Y[-(i + 1)] = yy2
            Z[-(i + 1)] = zz2
        elif d >= dr_left:
            # do some magic with reiterating the last step
            if flag_verbose:
                print( "too close. line step, dist", d, " dist left = ", dr_left)
            # xx, yy1 = xx_i, yy1_i
            # x2, y2 = x2_i, y2_i
            x = np.linspace(xx1, xx2, Nleft, endpoint=False, retstep=False)
            y = np.linspace(yy1, yy2, Nleft, endpoint=False, retstep=False)
            z = np.linspace(zz1, zz2, Nleft, endpoint=False, retstep=False)
            for j in range(Nleft):
                X[i+j+2], Y[i+j+2], Z[i+j+2] = x[j], y[j], z[j]
                if flag_verbose:
                    print( "steps ", x[j], y[j], z[j])
            break
        else:
            print( "error")
        if flag_verbose:
            print( "Nleft ", Nleft)
            print( "X  zero", N - np.count_nonzero(X))
            print( "Y  zero", N - np.count_nonzero(Y))
            print( "Z  zero", N - np.count_nonzero(Z))
    # if flag_verbose:
    #     print( X)
    #     print( Y)
    #     print( Z)
    # return X, Y, Z
    return np.column_stack((X[1:], Y[1:], Z[1:])), b

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # for xinit in [0, 2, 3, 5]:
    for xinit in [ 5]:
        for xend in [7]:
            x_1 = np.array([xinit, 1, 1])
            x_2 = np.array([xend, 1, 1])
            bounds = np.array([10,10,10])
            N = 25
            Coords,_ = chain_alng_line(x_1, x_2, N, bounds)
            # print( Coords)
            print( "chain given N ", N, " result size = ", len(Coords))

            # Coords_mn, _ = mons_alng_line(x_1, x_2, N, bounds)
            # print( "mons given N ", N, " result size = ", len(Coords_mn))

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

