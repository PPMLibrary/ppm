import sys
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy as np 

from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.patches import Polygon
from matplotlib.colors import *
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

progname = os.path.basename(sys.argv[0])

# simple color map
# should do it for now, but should be replaced by a generic implementation that
# picks arbitrary number of colors
cmap = {-1: 'k',\
        1 : 'r',\
        2 : 'b',\
        3 : 'g',\
        4 : 'y',\
        5 : 'c',\
        6 : 'm'}
    
def isreal(el):
    if el[-1] == -1: return False
    else: return True

def isghost(el):
    if el[-1] == -1: return True
    else: return False


class SubDomain(object):
    """Represents a subdomain of the domain decomposition
    
    This abstract class provides the basic functionality for storing subdomain
    information and computing the coordinates for the cuboid faces to be drawn
    by the matplotlib routines
    """

    def __init__(self):
        """Initialize a subdomain with default values."""
        self.minc = []
        self.maxc = []
        self.proc = 0
        self.bc = []
        self.halo = 0.0

    def __init__(self,minc,maxc,proc,bc,halo):
        """Initialize a subdomain with user provided information."""
        self.minc = minc[:]
        self.maxc = maxc[:]
        self.proc = proc
        self.bc = bc[:]
        self.halo = halo

    def __str__(self):
        return 'min: '+str(self.minc)+' max: '+str(self.maxc)
    
    def mkfaces(self,minc,maxc): pass
    
    def sub2faces(self): pass
    
    def gl2faces(self): pass


class SubDomain2(SubDomain):

    def __init__(self):
        SubDomain.__init__(self)
    
    def __init__(self,minc,maxc,proc,bc,halo):
        SubDomain.__init__(self,minc,maxc,proc,bc,halo)
   
    def sub2faces(self):
        """Return matplotlib compatible face coordinate information."""
        return self.mkfaces(self.minc,self.maxc)

    def mkfaces(self,minc,maxc):
        """Construct and return an array with the 4 corners of the rect."""
        x1 = minc[0]
        y1 = minc[1]
        x2 = maxc[0]
        y2 = maxc[1]
        return np.array([(x1,y1),(x2,y1),(x2,y2),(x1,y2)])

    def gl2faces(self):
        """Construct and return the rectangle including the ghostlayer."""
        gminc = [0.0]*2
        gmaxc = [0.0]*2
        gminc[0] = self.minc[0] - self.halo if self.bc[0] == 1 else self.minc[0]
        gmaxc[0] = self.maxc[0] + self.halo if self.bc[1] == 1 else self.maxc[0]
        gminc[1] = self.minc[1] - self.halo if self.bc[2] == 1 else self.minc[1]
        gmaxc[1] = self.maxc[1] + self.halo if self.bc[3] == 1 else self.maxc[1]
        return self.mkfaces(gminc,gmaxc)

class SubDomain3(SubDomain):
    
    def __init__(self):
        SubDomain.__init__(self)
    
    def __init__(self,minc,maxc,proc,bc,halo):
        SubDomain.__init__(self,minc,maxc,proc,bc,halo)

    def sub2faces(self):
        """Return matplotlib compatible 3D cuboid faces."""
        return self.mkfaces(self.minc,self.maxc)

    def mkfaces(self,minc,maxc):
        """Construct and return an array containing cuboid faces."""
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
    
    def gl2faces(self):
        """Construct and return the 3D cuboid including the ghostlayer."""
        gminc = [0.0]*3
        gmaxc = [0.0]*3
        gminc[0] = self.minc[0] - self.halo if self.bc[0] == 1 else self.minc[0]
        gmaxc[0] = self.maxc[0] + self.halo if self.bc[1] == 1 else self.maxc[0]
        gminc[1] = self.minc[1] - self.halo if self.bc[2] == 1 else self.minc[1]
        gmaxc[1] = self.maxc[1] + self.halo if self.bc[3] == 1 else self.maxc[1]
        gminc[2] = self.minc[2] - self.halo if self.bc[4] == 1 else self.minc[2]
        gmaxc[2] = self.maxc[2] + self.halo if self.bc[5] == 1 else self.maxc[2]
        return self.mkfaces(gminc,gmaxc)

class DomainDecomp(object):
    """Represents a 2D or 3D ppm domain decomposition.
    
    This data type provides functionality to represent the ppm domain
    decomposition. Further, it returns a list of faces or lines to be drawn by
    matplotlib routines."""
    
    def __init__(self,dim):
        """Initialize an empty domain decomposition and set dimensionality."""
        self.subs = []
        self.dim = dim

    def addSub(self,minc,maxc,proc,bc,halo):
        """Add a new subdomain."""
        if self.dim == 2:
            newSub = SubDomain2(minc,maxc,proc,bc,halo)
        elif self.dim == 3:
            newSub = SubDomain3(minc,maxc,proc,bc,halo)
        self.subs.append(newSub)

    def getFaces(self):
        """Get all subdomain faces."""
        faces = []
        for sub in self.subs:
            faces.append((sub.sub2faces(),sub.proc))
        return faces
    
    def getGhostFaces(self):
        """Get all ghostlayer faces."""
        faces = []
        for sub in self.subs:
            faces.append(sub.gl2faces())
        return faces



class Particles(object):
    """Represents a set of particles."""

    def __init__(self):
        self.xp = []
        self.c = []
        self.dim = 2

    def __init__(self,dim):
        """Initialize an empty set of particles in a dim-dimensional space."""
        self.xp = []
        self.c = []
        self.dim = dim
 
    def add(self,p,c):
        """Add a particle and its color tag."""
        self.xp.append(p)
        self.c.append(c)

    def x(self):
        """Return the x coordinate of all particles."""
        return [p[0] for p in self.xp]
    
    def y(self):
        """Return the y coordinate of all particles."""
        return [p[1] for p in self.xp]
    
    def z(self):
        """Return the z coordinate of all particles."""
        if self.dim == 3:
            return [p[2] for p in self.xp]
        else:
            return None


class AppForm(QMainWindow):
    """The main program and UI class."""
    
    def __init__(self, parent=None):
        self.halo = 0.0
        self.dim = 2
        self.decomp = None 
        self.particles = None
    
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('PPM Debug Utility')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        

    def load_data(self,subf,datf):
        """Read domain decomp and particles from file and store them."""
        l1 = subf.readline()
        self.dim = float(l1.strip())
        l2 = subf.readline()
        self.halo = float(l2.strip())
        self.decomp = DomainDecomp(self.dim)
    
        if self.dim == 2:
            for l in subf:
                r = l.strip().split()
                min_sub = [float(r[0]),float(r[1])]
                max_sub = [float(r[2]),float(r[3])]
                proc = int(r[4])
                bc = [int(r[5]),int(r[6]),int(r[7]),int(r[8])]
                self.decomp.addSub(min_sub,max_sub,proc,bc,self.halo)
        elif self.dim == 3:
            for l in subf:
                r = l.strip().split()
                min_sub = [float(r[0]),float(r[1]),float(r[2])]
                max_sub = [float(r[3]),float(r[4]),float(r[5])]
                proc = int(r[6])
                bc = [int(r[7]),int(r[8]),int(r[9]),\
                        int(r[10]),int(r[11]),int(r[12])]
                self.decomp.addSub(min_sub,max_sub,proc,bc,self.halo)
        
        if datf:
            self.particles = Particles(self.dim)
            if self.dim == 2:
                for l in datf:
                    r = l.strip().split()
                    xp = [float(r[0]),float(r[1])]
                    c = int(r[2])
                    self.particles.add(xp,c)
            elif self.dim == 3:
                for l in datf:
                    r = l.strip().split()
                    xp = [float(r[0]),float(r[1]),float(r[2])]
                    c = int(r[3])
                    self.particles.add(xp,c)
        else:
            self.particles = None
  
    def on_pick(self):
        """ pick event handler for mpl canvas."""
        pass
    
    def on_wheel(self,event):
        zoom_factor = 0.8
        if event.step > 0:
            zoom = zoom_factor/event.step
        else:
            zoom = (-1*event.step)/zoom_factor
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        xdist = xlim[1]-xlim[0]
        ydist = ylim[1]-ylim[0]
        zoom_xdist = (xdist*zoom)/2.0
        zoom_ydist = (ydist*zoom)/2.0

        if self.dim == 2:
            zoom_center = (event.xdata,event.ydata)
            new_xlim = (zoom_center[0] - zoom_xdist, \
                        zoom_center[0] + zoom_xdist)
            new_ylim = (zoom_center[1] - zoom_ydist, \
                        zoom_center[1] + zoom_ydist)
        elif self.dim == 3:
            zlim = self.axes.get_zlim3d()
            zdist = (zlim[1]-zlim[0])*zoom
            zoom_zdist = (zdist*zoom)/2.0
            zoom_center = [xlim[0]+xdist/2.0, \
                           ylim[0]+ydist/2.0, \
                           zlim[0]+zdist/2.0]
            new_xlim = (zoom_center[0] - zoom_xdist, \
                        zoom_center[0] + zoom_xdist)
            new_ylim = (zoom_center[1] - zoom_ydist, \
                        zoom_center[1] + zoom_ydist)
            new_zlim = (zoom_center[1] - zoom_zdist, \
                        zoom_center[1] + zoom_zdist)
            self.axes.set_zlim3d(new_zlim)
        
        self.axes.set_xlim(new_xlim)
        self.axes.set_ylim(new_ylim)
        
        self.canvas.draw()
        if self.dim == 3:
            self.axes.mouse_init()

    def on_draw(self):
        """ Redraws the figure."""
        self.axes.clear()
        if self.dim == 2:
            for el in self.decomp.getGhostFaces():
                self.plotgl2(el)
            for el,cpu in self.decomp.getFaces():
                self.plotsub2(el,cpu)
            self.plotdat2(self.particles.x(),\
                    self.particles.y(),\
                    self.particles.c)
        elif self.dim == 3:
            for el in self.decomp.getGhostFaces():
                self.plotgl3(el)
            for el,cpu in self.decomp.getFaces():
                self.plotsub3(el,cpu)
            self.plotdat3(self.particles.x(),\
                    self.particles.y(),\
                    self.particles.z(),\
                    self.particles.c)
        else: pass
        self.axes.set_xlabel('x')
        self.axes.set_ylabel('y')
        if self.dim == 3:
            self.axes.set_zlabel('z')
        
        self.canvas.draw()
        if self.dim == 3:
            self.axes.mouse_init()


    def open_data(self):
        """Open subdomain and data files chosen by the user."""
        file_choices = "*.sub"
        
        subfpath = unicode(QFileDialog.getOpenFileName(self, 
                        'Open file', '', 
                        file_choices))
        if not subfpath:
            return
        datfpath = subfpath[:-4]+'.dat'

        self.statusBar().showMessage('Opening %s' % subfpath, 2000)
        
        subf = None
        datf = None

        try:
            subf = open(subfpath,'r')
        except: pass
        try:
            datf = open(datfpath,'r')
        except: pass
        self.load_data(subf,datf)
        subf.close()
        datf.close()
        if self.dim == 2:
            self.axes = self.fig.add_subplot(111)
        else:
            self.axes = Axes3D(self.fig)
        self.on_draw()

    
    def save_plot(self):
        """Save the current plot to a png file."""
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
   

    def on_about(self):
        msg = """ PPM Debug Utility (c) 2011 MOSAIC Group
        Created by Omar Awile.
        
        """
        QMessageBox.about(self, "About ppmdbg", msg.strip())

    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        #self.axes = Axes3D(self.fig)
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('scroll_event', self.on_wheel)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # Other GUI controls
        # 
        self.draw_button = QPushButton("&Draw")
        self.connect(self.draw_button, SIGNAL('clicked()'), self.on_draw)
        
        #
        # Layout with box sizers
        # 
        hbox = QHBoxLayout()
        
        for w in [  self.draw_button]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignVCenter)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(hbox)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
       
        open_data_action = self.create_action("&Open data file",
                shortcut="Ctrl+O",slot=self.open_data,
                tip="Open set of debug files")
        save_plot_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (open_data_action,save_plot_action, None, quit_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')

    def create_status_bar(self):
        self.status_text = QLabel("status")
        self.statusBar().addWidget(self.status_text, 1)
    
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)
    
    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def plotsub3(self,f,cpu):
        """Plot a 3D subdomain."""
        nc = len(cmap.keys())
        for i in range(6):
            try:
                self.axes.plot_surface(f[i][0],f[i][1],f[i][2],alpha=0.05,\
                        color=cmap[cpu%(nc-1)+1]) 
            except KeyError:
                print "invalid color tag"
    
    def plotgl3(self,gl):
        """Plot a 3D subdomain ghostlayer."""
        for i in range(6):
            self.axes.plot_surface(gl[i][0],gl[i][1],gl[i][2],alpha=0.02,\
                    color='k',linewidth=0) 
    
    def plotsub2(self,f,cpu):
        """Plot a 2D subdomain."""
        nc = len(cmap.keys())
        try:
            p = Polygon(f,alpha=0.05,color=cmap[cpu%(nc-1)+1],linewidth=0)
        except KeyError:
            print "invalid color tag"
        self.axes.add_patch(p)
        p = Polygon(f,fill=False,linewidth=1,ec='k')
        self.axes.add_patch(p)
    
    def plotgl2(self,gl):
        """Plot a 2D subdomain ghostlayer."""
        p = Polygon(gl,alpha=0.01,color='k')
        self.axes.add_patch(p)
        p = Polygon(gl,fill=False,linewidth=0.4,linestyle='dashed',ec='k')
        self.axes.add_patch(p)
   
    
    def plotdat2(self,x,y,tag):
        """Plot 2D particle positions."""
        nc = len(cmap.keys())
        try:
            rx,ry,rtag = zip(*filter(isreal,zip(x,y,tag)))
            self.axes.scatter(rx,ry,s=5,c=[cmap[t%(nc-1)+1] for t in \
                rtag],linewidths=0)
        except KeyError:
            print "invalid color tag"
        try:
            gx,gy,gtag = zip(*filter(isghost,zip(x,y,tag)))
            self.axes.scatter(gx,gy,s=5,c=[cmap[t] for t in \
                gtag],linewidths=0,alpha=0.6,zorder=10)
        except KeyError:
            print "invalid color tag"
        except ValueError:
            print "no ghosts"

    
    def plotdat3(self,x,y,z,tag):
        """Plot 3D particle positions."""
        nc = len(cmap.keys())
        try:
            rx,ry,rz,rtag = zip(*filter(isreal,zip(x,y,z,tag)))
            self.axes.scatter(rx,ry,rz,s=10,c=[cmap[t%(nc-1)+1] for t in \
                rtag],linewidths=0)
        except KeyError:
            print "invalid color tag"
        try:
            gx,gy,gz,gtag = zip(*filter(isghost,zip(x,y,z,tag)))
            self.axes.scatter(gx,gy,gz,s=10,c=[cmap[t] for t in \
                gtag],linewidths=0,alpha=0.6)
        except KeyError:
            print "invalid color tag"
        except ValueError:
            print "no ghosts"


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
