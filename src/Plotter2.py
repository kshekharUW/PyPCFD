'''
Created on Jun 13, 2018
Modified on June 24, 2018 for separate images

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
import matplotlib.cm as cm
import os


class Plotter(object):
    '''
    variables:
        self.particlesPresent     # boolean
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
        def setData(self, nodes)
        def setParticleData(self, particles)
        def setCellFluxData(self, cells)
        def plotCellFlux(self, time)
        
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.IMAGE_COUNTER = -1
        self.particlesPresent = False
        
        self.width  = 1.0
        self.height = 1.0
        
        self.setGrid(self.width, self.height, 10, 10)

        if not os.path.isdir("images"):
            os.mkdir("images")


    def refresh(self, time=-1):
        
        self.IMAGE_COUNTER += 1
        
        speed = np.sqrt(self.Vx*self.Vx + self.Vy*self.Vy )
        
        fig = plt.figure(figsize=(10, 9))
        #gs = gridspec.GridSpec(nrows=1, ncols=1, height_ratios=[1, 1], width_ratios=[1,1.3])
        
        #  pressure field
        #ax0 = fig.add_subplot(gs[1, 1])
        ax0 = fig.gca()
        try:
            contour = ax0.contourf(self.X, self.Y, self.P, cmap='autumn')
            fig.colorbar(contour, ax=ax0)
            if (time>=0.0):
                ax0.set_title('Pressure at t={:08.5f}s'.format(time))
            else:
                ax0.set_title('Pressure')
        except:
            ax0.set_title('Pressure Failed at t={:08.5f}s'.format(time))
        
        imageName = "Pressure{:04d}.png".format(self.IMAGE_COUNTER)
        plt.savefig("images/"+imageName)
        
        plt.clf()
        
        # Varying color along a streamline
        #ax1 = fig.add_subplot(gs[1, 0])
        ax1 = fig.gca()
        try:
            force = ax1.quiver(self.X, self.Y, self.Fx, self.Fy, cmap='autumn')
            ax1.quiverkey(force, 0.9*self.width, 1.05*self.height, 1.0, r'$1.0\,N$',labelpos='E',coordinates='axes')
            #fig.colorbar(force, ax=ax1)
            if (time>=0.0):
                ax1.set_title('Nodal Forces at t={:08.5f}s'.format(time))
            else:
                ax1.set_title('Nodal Forces')
        except:
            ax1.set_title('Nodal Forces Failed at t={:08.5f}s'.format(time))
        
        imageName = "Forces{:04d}.png".format(self.IMAGE_COUNTER)
        plt.savefig("images/"+imageName)
        
        plt.clf()
        
        
        #  Varying line width along a streamline
        #ax2 = fig.add_subplot(gs[0, 0])
        ax2 = fig.gca()
        try:
            vecs = ax2.quiver(self.X, self.Y, self.Vx, self.Vy, cmap='autumn')
            ax2.quiverkey(vecs, 0.9*self.width, 1.05*self.height, 1.0, r'$1.0\,\frac{m}{s}$',labelpos='E',coordinates='axes')
            #fig.colorbar(vecs, ax=ax2)
            if (time>=0.0):
                ax2.set_title('velocity at t={:08.5f}s'.format(time))
            else:
                ax2.set_title('velocity')
        except:
            ax2.set_title('velocity Failed at t={:08.5f}s'.format(time))
        
        imageName = "Velocity{:04d}.png".format(self.IMAGE_COUNTER)
        plt.savefig("images/"+imageName)
        
        plt.clf()
        
        
        #ax3 = fig.add_subplot(gs[0, 1])
        ax3 = fig.gca()
        try:
            seed_points = np.array([ self.tracerPoints[0].flatten(), self.tracerPoints[1].flatten() ])
            
            #vecs  = ax3.quiver(    self.X, self.Y, self.Vx, self.Vy, cmap='autumn')
            strm3 = ax3.streamplot(self.X, self.Y, self.Vx, self.Vy, 
                                   color=speed, linewidth=1, cmap='autumn',
                                   density=2, integration_direction='forward',
                                   start_points=seed_points.T)
            fig.colorbar(strm3.lines)
            if (time>=0.0):
                ax3.set_title('Streamlines at t={:08.5f}s'.format(time))
            else:
                ax3.set_title('Streamlines')
        except:
            ax3.set_title('Streamlines Failed at t={:08.5f}s'.format(time))
        
        # Displaying the starting points with blue symbols.
        ###ax3.plot(self.tracerPoints[0], self.tracerPoints[1], 'bo', markersize=2)
        ax3.axis((0.0, self.width, 0.0, self.height))
        
        imageName = "Stream{:04d}.png".format(self.IMAGE_COUNTER)
        plt.savefig("images/"+imageName)
        
        plt.clf()
        
        
        if (self.particlesPresent):
            #  Varying line width along a streamline
            #ax2 = fig.add_subplot(gs[0, 0])
            ax4 = fig.gca()
            try:
                vecs = ax4.quiver(self.ParticleX, self.ParticleY, self.ParticleVx, self.ParticleVy, cmap='autumn')
                ax4.quiverkey(vecs, 0.9*self.width, 1.05*self.height, 1.0, r'$1.0\,\frac{m}{s}$',labelpos='E',coordinates='axes')
                
                points = ax4.plot(self.ParticleX, self.ParticleY, 'bo', markersize=2)
                #fig.colorbar(vecs, ax=ax2)
                if (time>=0.0):
                    ax4.set_title('particle velocity at t={:08.5f}s'.format(time))
                else:
                    ax4.set_title('particle velocity')
            except:
                ax4.set_title('particle velocity Failed at t={:08.5f}s'.format(time))
                
            #ax4.axis((0, 1, 0, 1))
            
            imageName = "ParticleVelocity{:04d}.png".format(self.IMAGE_COUNTER)
            plt.savefig("images/"+imageName)

            plt.clf()


        
        plt.close()
        self.plotCellFlux(time)
        
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
        
        xp = self.width*(0.5+np.arange(nCellsX)) / nCellsX
        yp = self.height*(0.5+np.arange(nCellsY))/ nCellsY
        refPtsX, refPtsY = np.meshgrid(xp,yp)
        
        self.tracerPoints = [refPtsX,refPtsY]
        
    def setData(self, nodes):
        self.particlesPresent = False
        
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
        
    def setParticleData(self, particles):
        if (len(particles) == 0):
            return
        
        self.particlesPresent = True
        
        self.ParticleX  = []
        self.ParticleY  = []
        self.ParticleVx = []
        self.ParticleVy = []
        
        for p in particles:
            pos = p.position()
            vel = p.velocity()
            self.ParticleX.append(pos[0])
            self.ParticleY.append(pos[1])
            self.ParticleVx.append(vel[0])
            self.ParticleVy.append(vel[1])

    def setCellFluxData(self, cells):
        ncellsX = self.nNodesX-1
        ncellsY = self.nNodesY-1
        self.cellFlux = np.zeros((ncellsX, ncellsY))
        for theCell in cells:
            cellGridCoords = theCell.getCellGridCoordinates()
            self.cellFlux[cellGridCoords[0]][cellGridCoords[1]] = theCell.getFlux()


    def plotCellFlux(self, time):
        fig, ax = plt.subplots()
        img = ax.imshow(np.flipud(self.cellFlux), cmap=cm.jet, interpolation='nearest', vmin=-0.5, vmax=0.5)
        # fig.colorbar(img)
        # Loop over data dimensions and create text annotations.
        for i in range(self.nNodesX-1):
            for j in range(self.nNodesY-1):
                strToDisplay = "{:10.4f}".format(self.cellFlux[i, j])
                text = ax.text(j, i, strToDisplay,
                               ha="center", va="center", color="k")

        if (time >= 0.0):
            ax.set_title('cell flux at t={:08.5f}s'.format(time))
        else:
            ax.set_title('cell flux')

        ax.set_xticks([])
        ax.set_yticks([])
        imageName = "CellFlux{:04d}.png".format(self.IMAGE_COUNTER)
        plt.savefig("images/" + imageName)

        plt.close()
