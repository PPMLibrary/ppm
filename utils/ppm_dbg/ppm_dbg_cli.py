from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Polygon
from matplotlib.colors import *
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import sys

# should do it for now, but should be replaced by a generic implementation that
# picks arbitrary number of colors
cmap = {-1: 'k',\
        1 : 'r',\
        2 : 'b',\
        3 : 'g',\
        4 : 'y',\
        5 : 'c',\
        6 : 'm'}

halo = 0.0

def sub2rect(minc,maxc):
    x1 = minc[0]
    y1 = minc[1]
    x2 = maxc[0]
    y2 = maxc[1]
    return np.array([(x1,y1),(x2,y1),(x2,y2),(x1,y2)])

def gl2rect(minc,maxc,bc,halo):
    gminc = [0.0]*2
    gmaxc = [0.0]*2
    gminc[0] = minc[0] - halo if bc[0] == 1 else minc[0]
    gmaxc[0] = maxc[0] + halo if bc[1] == 1 else maxc[0]
    gminc[1] = minc[1] - halo if bc[2] == 1 else minc[1]
    gmaxc[1] = maxc[1] + halo if bc[3] == 1 else maxc[1]
    return sub2rect(gminc,gmaxc)


def sub2cube(minc,maxc):
    x1 = minc[0]
    y1 = minc[1]
    z1 = minc[2]
    x2 = maxc[0]
    y2 = maxc[1]
    z2 = maxc[2]

    # front face
    f1x = np.array([x1,x2])
    f1y = np.array([[y1,y1],[y1,y1]])
    f1z = np.array([z1,z2])
    f1x,f1z = np.meshgrid(f1x,f1z)

    # back face
    f2x = np.array([x1,x2])
    f2y = np.array([[y2,y2],[y2,y2]])
    f2z = np.array([z1,z2])
    f2x,f2z = np.meshgrid(f2x,f2z)

    # left face
    f3x = np.array([[x1,x1],[x1,x1]])
    f3y = np.array([y1,y2])
    f3z = np.array([z1,z2])
    f3y,f3z = np.meshgrid(f3y,f3z)


    # right face
    f4x = np.array([[x2,x2],[x2,x2]])
    f4y = np.array([y1,y2])
    f4z = np.array([z1,z2])
    f4y,f4z = np.meshgrid(f4y,f4z)

    # bottom face
    f5x = np.array([x1,x2])
    f5y = np.array([y1,y2])
    f5z = np.array([[z1,z1],[z1,z1]])
    f5x,f5y = np.meshgrid(f5x,f5y)

    # top face
    f6x = np.array([x1,x2])
    f6y = np.array([y1,y2])
    f6z = np.array([[z2,z2],[z2,z2]])
    f6x,f6y = np.meshgrid(f6x,f6y)

    return [(f1x,f1y,f1z),(f2x,f2y,f2z),(f3x,f3y,f3z),(f4x,f4y,f4z),\
            (f5x,f5y,f5z),(f6x,f6y,f6z)]

def gl2cube(minc,maxc,bc,halo):
    gminc = [0.0]*3
    gmaxc = [0.0]*3
    gminc[0] = minc[0] - halo if bc[0] == 1 else minc[0]
    gmaxc[0] = maxc[0] + halo if bc[1] == 1 else maxc[0]
    gminc[1] = minc[1] - halo if bc[2] == 1 else minc[1]
    gmaxc[1] = maxc[1] + halo if bc[3] == 1 else maxc[1]
    gminc[2] = minc[2] - halo if bc[4] == 1 else minc[2]
    gmaxc[2] = maxc[2] + halo if bc[5] == 1 else maxc[2]
    return sub2cube(gminc,gmaxc)

def plotsub3(ax,f,cpu):
    nc = len(cmap.keys())
    for i in range(6):
        ax.plot_surface(f[i][0],f[i][1],f[i][2],alpha=0.05,color=cmap[cpu%(nc+1)+1])

def plotgl3(ax,gl):
    for i in range(6):
        ax.plot_surface(gl[i][0],gl[i][1],gl[i][2],alpha=0.02,\
                color='k',linewidth=0)

def plotsub2(ax,f,cpu):
    nc = len(cmap.keys())
    p = Polygon(f,alpha=0.05,color=cmap[cpu%(nc+1)+1],linewidth=0)
    ax.add_patch(p)
    p = Polygon(f,fill=False,linewidth=1,ec='k')
    ax.add_patch(p)

def plotgl2(ax,gl):
    p = Polygon(gl,alpha=0.01,color='k')
    ax.add_patch(p)
    p = Polygon(gl,fill=False,linewidth=0.4,linestyle='dashed',ec='k')
    ax.add_patch(p)

def plotdat2(ax,x,y,tag):
    ax.scatter(x,y,s=5,c=[cmap[t] for t in tag],linewidths=0)

def plotdat3(ax,x,y,z,tag):
    ax.scatter(x,y,z,s=10,c=[cmap[t] for t in tag],linewidths=0)

def main():
    if len(sys.argv) > 1:
      subfilen = sys.argv[1]
      datfilen = sys.argv[2]

      fig = plt.figure()
      subfile = open(subfilen)
      l1 = subfile.readline()
      dim = float(l1.strip())
      print dim
      l2 = subfile.readline()
      halo = float(l2.strip())

      if dim == 2:
        ax = fig.add_subplot(111)
        for l in subfile:
          r = l.strip().split()
          min_sub = [float(r[0]),float(r[1])]
          max_sub = [float(r[2]),float(r[3])]
          proc = int(r[4])
          bc = [int(r[5]),int(r[6]),int(r[7]),int(r[8])]
          if (halo > 0.0):
            glfaces = gl2rect(min_sub,max_sub,bc,halo)
            plotgl2(ax,glfaces)
            faces = sub2rect(min_sub,max_sub)
            plotsub2(ax,faces,proc)
      elif dim == 3:
        ax = Axes3D(fig)
        for l in subfile:
          r = l.strip().split()
          min_sub = [float(r[0]),float(r[1]),float(r[2])]
          max_sub = [float(r[3]),float(r[4]),float(r[5])]
          proc = int(r[6])
          bc = [int(r[7]),int(r[8]),int(r[9]),\
            int(r[10]),int(r[11]),int(r[12])]
          if (halo > 0.0):
            glfaces = gl2cube(min_sub,max_sub,bc,halo)
            plotgl3(ax,glfaces)
            faces = sub2cube(min_sub,max_sub)
            plotsub3(ax,faces,proc)

      subfile.close()
      datfile = open(datfilen)
      x = []
      y = []
      if dim == 3:
        z = []

      c = []

      if dim == 2:
        for l in datfile:
          r = l.strip().split()
          x.append(float(r[0]))
          y.append(float(r[1]))
          c.append(int(r[2]))
          plotdat2(ax,x,y,c)
      elif dim == 3:
        for l in datfile:
          r = l.strip().split()
          x.append(float(r[0]))
          y.append(float(r[1]))
          z.append(float(r[2]))
          c.append(int(r[3]))
          plotdat3(ax,x,y,z,c)

      datfile.close()

      ax.set_xlabel('x')
      ax.set_ylabel('y')
      if dim == 3:
        ax.set_zlabel('z')

      plt.show()

if __name__ == "__main__":
    main()
