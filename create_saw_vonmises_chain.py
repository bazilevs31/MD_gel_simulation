"""Summary
"""
import numpy as np
import distances

def cart2spher(x_):
    """return spheraical coordinates of a cartesian vector

    Args:
        x_ (TYPE): Description
    """
    x, y, z = x_[0], x_[1], x_[2]
    rho = np.sqrt(x**2 + y**2 + z**2)
    theta = z/float(rho)
    phi = np.arctan2(y, x)
    return(rho, phi, theta)

def dist(v1, v2):
    """returns the dist between two vectors

    Args:
        v1 (TYPE): Description
        v2 (TYPE): Description
    """
    # return np.linalg.norm(x-y)
    a = v1 - v2
    return np.sqrt(np.einsum('i,i', a, a))

def get_kappa(x, y, nrem, N, b=1.):
    """get an emperic kappa parameter for chain creation
    # Note higher "kappas" (second arg) result in a _narrower_ distribution

    Args:
        x (TYPE): Description
        y (TYPE): Description
        nrem (TYPE): Description
        N (TYPE): Description
        b (float, optional): Description
    """
    d = dist(x, y)
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

def get_vect(x, k, phi, theta, b=1.):
    """create a random vector based on a certain direction
    from a given point
    uses a von mises distribution

    Args:
        x (TYPE): Description
        k (TYPE): Description
        phi (TYPE): Description
        theta (TYPE): Description
        b (float, optional): Description
    """
    x_new = np.zeros(3)

    alpha = np.random.vonmises(phi, k)
    theta = np.random.vonmises(theta, k)
    # print "alpha, theta ", alpha, theta

    x_new[0] = x[0] + b * np.cos(alpha) * np.cos(theta)
    x_new[1] = x[1] + b * np.sin(alpha) * np.cos(theta)
    x_new[2] = x[2] + b * np.sin(theta)
    return x_new





def create_chain(x1, x2, N, bounds, b=1.):
    """create a random chain using a self-drunk walk (von mises distrib)
    with a variable width

    Args:
        x1 (TYPE): Description
        x2 (TYPE): Description
        N (TYPE): Description
        bounds (TYPE): Description
        b (float, optional): Description
    """
    x1_i = np.zeros_like(x1)
    x2 = distances.get_pbc(x1, x2, bounds)
    Coord = np.zeros((N, 3))
    Coord[0,:] = x1[:]
    Coord[-1,:] = x2[:]
    for nsteps in range(1):
        for i in range(1, N-1):
            _, phi_direction, theta_direction = cart2spher(x2-x1)
            nrem = N - (i+2)
            kappa = get_kappa(x1, x2 , nrem, N)
            x1_i = get_vect(x1, kappa, phi_direction, theta_direction, b)
            # print "nrem = ", nrem, "dist = ", dist(x1, x1_i)
            # print kappa, x1_i, x1

            x1 = np.copy(x1_i)
            Coord[i,:] = x1[:]
        # print "nrem = ", nrem, "dist = ", dist(x1, x2)
    return Coord[1:-1,:], b

def writexyz(X, Y, Z, filename='mol_default.xyz'):
        """writes an xyz file that could vizualized using vmd

        Args:
            X (TYPE): Description
            Y (TYPE): Description
            Z (TYPE): Description
            filename (str, optional): Description

        Returns:
            TYPE: Description
        """
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

if __name__ == '__main__':
    N = 100
    b = 1.

    x1 = np.array([0,0,0])
    x2 = np.array([15,0,0])
    bounds = np.array([10, 10, 10])
    Coord = create_chain(x1, x2, N, bounds)
    X, Y, Z = Coord[:,0], Coord[:,1], Coord[:,2]
    writexyz(X, Y, Z)
