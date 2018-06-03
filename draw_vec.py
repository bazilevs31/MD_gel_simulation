import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

# #draw cube
# r = [-1, 1]
# for s, e in combinations(np.array(list(product(r,r,r))), 2):
#     if np.sum(np.abs(s-e)) == r[1]-r[0]:
#         ax.plot3D(*zip(s,e), color="b")

# #draw sphere
# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# x=np.cos(u)*np.sin(v)
# y=np.sin(u)*np.sin(v)
# z=np.cos(v)
# ax.plot_wireframe(x, y, z, color="r")

#draw a point
ax.scatter([0],[0],[0],color="g",s=100)

r0 = np.array([1.,1.,1.])
r1 = np.array([2.,2.,1.])
alpha = np.pi/6.
rtr = r0+ (r1-r0)/(2*np.cos(alpha))

Coord =np.vstack((r0, r1, rtr))

ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)

for i in range(len(Coord[:,0])):
    a = Arrow3D([0,Coord[i,0]],[0,Coord[i,1]],[0,Coord[i,2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
    ax.add_artist(a)
plt.show()
