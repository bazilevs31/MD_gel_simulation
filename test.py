import numpy as np
import numba

import time



# @profile
@numba.jit(nopython=True)
def get_cryst(bonds, alignment_vector, threshold):
    """
    input: bonds, alignment_vector, threshold
    output: cryst


    bonds : numpy array of (N,3) of normalized bonds vectors of a selected
    domain
    alignment_vector: eigen vector, which corresponds to orientation of
    selected bonds, by which they will be compared with
    cryst: crystallinity of the region, analyzed by alignment

    then calculate the number of bonds that are aligned
    to the average direction alignment_vector
    then look at the percentage of the bonds which are within certain threshold
    """
    Nbonds = bonds.shape[0]
    k = 0
    cryst, P2, cos = 0.0, 0.0, 0.0
    for i in range(Nbonds):
        cos = 0.
        tmp = 0.0
        for j in range(3):
            tmp = alignment_vector[j]*bonds[i, j]
            cos += tmp
        # cos = np.dot(alignment_vector, bonds[i])
        # print cos
        # cos = 0.
        P2 = (3.*cos**2.0 - 1.)/2.
        if P2 > threshold:
            k += 1
            # print "k is larger", k
    # print "k = ", k
    if k > 0:
        # print "k = ", k
        cryst = float(k)/float(Nbonds)
        # print "crystallinity" , cryst
    elif k == 0:
        cryst = 0.
    return cryst


N = 100000
bonds = np.zeros((N, 3))
for i in range(N):
    bonds[i, 0] = 1.
bonds += 0.15*np.random.rand(N, 3)
# bonds = np.eye(3)
# bonds = np.eye(3)
normb = np.linalg.norm(bonds, axis=1)
bonds /= normb[:, None]

# norm = np.linalg.norm(bonds)
# bonds /= norm
print bonds

# bonds = np.ones((1000, 3))
# alignment_vector = np.te, axis=1)
alignment_vector = np.array([1.0,0.,0.])
norm = np.linalg.norm(alignment_vector)
# alignment_vector /= norm[:, None]
alignment_vector /= norm
print alignment_vector
t0 = time.time()

result = get_cryst(bonds, alignment_vector, 0.98)
# get_cryst(bonds, alignment_vector)
t1 = time.time()

print "total time", t1-t0
print "overall crystallinity", result
