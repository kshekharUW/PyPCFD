'''
Created on Jun 13, 2018

@author: pmackenz
'''
import numpy as np
from sys import platform

import matplotlib
if "win" in platform.lower():
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from matplotlib import cm


class Plotter(object):
    '''
    variables:
        self.height
        self.width
        self.nNodesX = nCellsX+1
        self.nNodesY = nCellsY+1
        self.Y
        self.X
        self.Vx
        self.Vy
        self.speed
        self.fig = plt.figure(figsize=(9, 9))
        self.gs
        self.tracerPoints = [[],[]]
    
    methods:
        def __init__(self)
        def safePlot(self, filename)
        def refresh(self, time=-1)
        def setGrid(self, nodes)
        
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.IMAGE_COUNTER = -1
        
        self.width  = 1.0
        self.height = 1.0
        
        self.setGrid(self.width, self.height, 10, 10)
        
    def safePlot(self, filename):
        self.fig.saveFig(filename)
        
    def refresh(self, time=-1):
        
        #self.Vx =  self.Y   # temporary solution for debugging only
        #self.Vy = -self.X   # temporary solution for debugging only
        
        speed = np.sqrt(self.Vx*self.Vx + self.Vy*self.Vy )
        
        fig = plt.figure(figsize=(10, 9))
        gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1, 1], width_ratios=[1,1.3])
        
        
        #  pressure field
        ax0 = fig.add_subplot(gs[1, 1])
        try:
            contour = ax0.contourf(self.X, self.Y, self.P, cmap='autumn')
            fig.colorbar(contour, ax=ax0)
            ax0.set_title('Pressure')
        except:
            ax0.set_title('Pressure Failed')
        
        # Varying color along a streamline
        ax1 = fig.add_subplot(gs[1, 0])
        try:
            force = ax1.quiver(self.X, self.Y, self.Fx, self.Fy, cmap='autumn')
            ax1.quiverkey(force, 0.9*self.width, 1.05*self.height, 1.0, r'$1.0\,N$',labelpos='E',coordinates='axes')
            #fig.colorbar(force, ax=ax1)
            ax1.set_title('Nodal Forces')
        except:
            ax1.set_title('Nodal Forces Failed')
        
        
        #  Varying line width along a streamline
        ax2 = fig.add_subplot(gs[0, 0])
        try:
            vecs = ax2.quiver(self.X, self.Y, self.Vx, self.Vy, cmap='autumn')
            ax2.quiverkey(vecs, 0.9*self.width, 1.05*self.height, 1.0, r'$1.0\,\frac{m}{s}$',labelpos='E',coordinates='axes')
            #fig.colorbar(vecs, ax=ax2)
            ax2.set_title('velocity')
        except:
            ax2.set_title('velocity Failed')
        
        
        ax3 = fig.add_subplot(gs[0, 1])
        try:
            seed_points = np.array([ self.tracerPoints[0].flatten(), self.tracerPoints[1].flatten() ])
            
            #vecs  = ax3.quiver(    self.X, self.Y, self.Vx, self.Vy, cmap='autumn')
            strm3 = ax3.streamplot(self.X, self.Y, self.Vx, self.Vy, 
                                   color=speed, linewidth=1, cmap='autumn',
                                   density=2, integration_direction='forward',
                                   start_points=seed_points.T)
            fig.colorbar(strm3.lines)
            ax3.set_title('Particle Trajectory')
        except:
            ax3.set_title('Particle Trajectory Failed')
        
        # Displaying the starting points with blue symbols.
        ###ax3.plot(self.tracerPoints[0], self.tracerPoints[1], 'bo', markersize=2)
        ax3.axis((0.0, self.width, 0.0, self.height))
        
        plt.tight_layout()
        self.IMAGE_COUNTER += 1
        imageName = "Stream{:03d}.png".format(self.IMAGE_COUNTER)
        plt.savefig(imageName)
        plt.close()
        
    def setGrid(self, width, height, nCellsX, nCellsY):
        self.height = height
        self.width  = width
        self.nNodesX = nCellsX+1
        self.nNodesY = nCellsY+1
        
        x = np.linspace(0,width ,(nCellsX+1))
        y = np.linspace(0,height,(nCellsY+1))
        
        self.X, self.Y = np.meshgrid(x, y)
        self.Vx = np.zeros_like(self.X)
        self.Vy = np.zeros_like(self.X)
        
        # define tracer points
        
        #xp = 0.5*self.width*(0.5+np.arange(2*nCellsX)) / nCellsX
        #yp = 0.5*self.height*(0.5+np.arange(2*nCellsY))/ nCellsY
        xp = self.width*(0.5+np.arange(nCellsX)) / nCellsX
        yp = self.height*(0.5+np.arange(nCellsY))/ nCellsY
        refPtsX, refPtsY = np.meshgrid(xp,yp)
        
        self.tracerPoints = [refPtsX,refPtsY]
        
    def setData(self, nodes):
        self.P  = np.zeros_like(self.X)
        self.Vx = np.zeros_like(self.X)
        self.Vy = np.zeros_like(self.X)
        self.Fx = np.zeros_like(self.X)
        self.Fy = np.zeros_like(self.X)
        
        for i in range(self.nNodesX):
            for j in range(self.nNodesY):
                node = nodes[i][j]
                self.P[j,i] = node.getPressure()
                vel = node.getVelocity()
                self.Vx[j,i] = vel[0]
                self.Vy[j,i] = vel[1]
                force = node.getForce()
                self.Fx[j,i] = force[0]
                self.Fy[j,i] = force[1]
        
        self.speed = np.sqrt(self.Vx*self.Vx + self.Vy*self.Vy)
        
        