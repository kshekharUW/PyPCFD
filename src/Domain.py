'''
Created on June 10, 2018

@author: pmackenz
''' 

from Node import *
from Cell import *
from Particle import *
from Plotter2 import *

from numpy import array, linspace, dot, tensordot, zeros, ones, outer, linspace, meshgrid, abs,\
    ceil
from numpy.linalg import solve
from xml.dom.minicompat import NodeList

class Domain(object):
    '''
    variable:
        self.width 
        self.height
        self.nCellsX 
        self.nCellsY 
        self.X 
        self.Y 
        self.hx        # cell size in x-direction
        self.hy        # cell size in y-direction
        self.nodes 
        self.cells 
        self.particles
        self.analysisControl
        self.Re
        self.v0
        self.time = 0.0
    
    methods:
        def __init__(self, width=1., height=1., nCellsX=2, nCellsY=2)
        def __str__(self)
        def setBoundaryConditions(self)
        def setAnalysis(self, doInit, solveVstar, solveP, solveVtilde, solveVenhanced, updatePosition, updateStress)
        def getAnalysisControl(self)
        def setStateN(self)
        def setParameters(self, Re, density, velocity)
        def runAnalysis(self, maxtime=1.0)
        def runSingleStep(self, dt=1.0)
        def initStep(self)
        def solveVstar(self, dt)
        def solveP(self, dt)
        def solveVtilde(self, dt)
        def solveVenhanced(self, dt)
        def updateParticleStress(self)
        def updateParticleMotion(self)
        def findCell(self, x)
        def createParticles(self, n, m)
    '''

    def __init__(self, width=1., height=1., nCellsX=2, nCellsY=2):
        '''
        Constructor
        '''
        self.width   = width
        self.height  = height
        self.nCellsX = nCellsX
        self.nCellsY = nCellsY
        
        self.hx = width/nCellsX        # cell size in x-direction
        self.hy = height/nCellsY       # cell size in y-direction
        
        self.time = 0.0
        
        #self.X = outer(ones(nCellsY+1), linspace(0.0, width, nCellsX+1))
        #self.Y = outer(linspace(0.0, height, nCellsY+1), ones(nCellsX+1))
        
        x = linspace(0,width ,(nCellsX+1))
        y = linspace(0,height,(nCellsY+1))
        
        self.X, self.Y = meshgrid(x, y, indexing='xy')
        
        self.Re  = 1.0
        self.rho = 1.0
        self.v0  = 0.0
        
        self.nodes = [ [ None for j in range(self.nCellsY+1) ] for i in range(self.nCellsX+1) ]
        id = -1
        
        for i in range(nCellsX+1):
            for j in range(nCellsY+1):
                id += 1
                theNode = Node(id,x[i],y[j])
                theNode.setGridCoordinates(i,j)
                self.nodes[i][j] = theNode
                
        self.cells = []
        id = -1
        hx = width / nCellsX
        hy = height / nCellsY
        
        for i in range(nCellsX):
            for j in range(nCellsY):
                id += 1
                newCell = Cell(id, hx, hy)
                theNodes = []
                theNodes.append(self.nodes[i][j])
                theNodes.append(self.nodes[i+1][j])
                theNodes.append(self.nodes[i+1][j+1])
                theNodes.append(self.nodes[i][j+1])
                newCell.SetNodes(theNodes)
                self.cells.append(newCell)     
        
        self.setParameters(self.Re, self.rho, self.v0)
        self.setAnalysis(False, True, True, True, True, True, False, False)
    
        self.plot = Plotter()
        self.plot.setGrid(width, height, nCellsX, nCellsY)
        
    def setBoundaryConditions(self):
        
        nCellsX = self.nCellsX
        nCellsY = self.nCellsY
        
        # define fixities
        for i in range(nCellsX+1):
            self.nodes[i][0].fixDOF(1, 0.0)
            self.nodes[i][nCellsY].fixDOF(1, 0.0)
            
            #self.nodes[i][0].fixDOF(0, 0.0)             # fully fixed
            #self.nodes[i][nCellsY].fixDOF(0, 0.0)       # fully fixed
            
            if (i>0 and i< nCellsX+1):
                self.nodes[i][nCellsY].fixDOF(0, self.v0)
        for j in range(nCellsY+1):
            self.nodes[0][j].fixDOF(0, 0.0)
            self.nodes[nCellsX][j].fixDOF(0, 0.0)
            
            #self.nodes[0][j].fixDOF(1, 0.0)             # fully xixed
            #self.nodes[nCellsX][j].fixDOF(1, 0.0)       # fully fixed
            
        
        
    def __str__(self):
        s = "==== D O M A I N ====\n"
        s += "Nodes:\n"
        for i in range(self.nCellsX+1):
            for j in range(self.nCellsY+1):
                s += str(self.nodes[i][j]) + "\n"
        s += "nCells:\n"
        for cell in self.cells:
            s += str(cell) + "\n"
        return s
        
    def setAnalysis(self, doInit, solveVstar, solveP, solveVtilde, solveVenhanced, updatePosition, updateStress, addTransient):
        self.analysisControl = {
            'doInit':doInit,
            'solveVstar':solveVstar,
            'solveP':solveP,
            'solveVtilde':solveVtilde,
            'solveVenhanced':solveVenhanced,
            'updatePosition':updatePosition,
            'updateStress':updateStress,
            'addTransient':addTransient
            }
        if (doInit and updatePosition and addTransient):
            print("INCONSISTENCY WARNING: transient active with updatePosition && doInit ")
        
    def getAnalysisControl(self):
        return self.analysisControl
    
    def setParameters(self, Re, density, velocity):
        
        if (self.hx < self.hy ):
            L = self.hx
        else:
            L = self.hy
            
        viscosity = density * velocity * L / Re
        
        self.Re = Re
        self.rho = density
        self.v0  = velocity
        self.mu  = viscosity
            
        self.setBoundaryConditions()
        
        for cell in self.cells:
            cell.setParameters(density, viscosity)
       
    def setInitialState(self):
        for nodeList in self.nodes:
            for node in nodeList:
                node.wipe()
        
        for cell in self.cells:
            cell.mapMassToNodes()
        
        # initial condition at nodes define v*, not v
        for i in range(self.nCellsX+1):
            for j in range(self.nCellsY):
                self.nodes[i][j].setVelocity(zeros(2))
            self.nodes[i][self.nCellsY].setVelocity(array([self.v0 ,0.0]))
        # fix the top corner nodes
        self.nodes[0][self.nCellsY].setVelocity(zeros(2))
        self.nodes[self.nCellsX][self.nCellsY].setVelocity(zeros(2))
        
        # now find pressure for a fictitious time step dt = 1.0
        self.solveP(1.0)
        
        # compute \tilde v
        self.solveVtilde(1.0)
        
        # initial conditions are now set
        self.plot.setData(self.nodes)
        self.plot.refresh(self.time)
     
    def runAnalysis(self, maxtime=1.0):
        
        # find ideal timestep using CFL=0.5
        dt = self.getTimeStep(1.0)
        if (dt > (maxtime - self.time)):
            dt = (maxtime - self.time)
        if (dt < (maxtime - self.time)):
            nsteps = ceil((maxtime - self.time)/dt)
            if (nsteps>50):
                nsteps= 50
            dt = (maxtime - self.time) / nsteps
        
        #
        while (self.time < maxtime-0.1*dt):
            self.runSingleStep(self.time, dt)
            self.time += dt 
            
        self.plot.setData(self.nodes)
        self.plot.refresh(self.time)
        
    
    def runSingleStep(self, time=0.0, dt=1.0):
        
        print("starting at t_n = {:.3f}, time step \u0394t = {}, ending at t_(n+1) = {:.3f}".format(time, dt, time+dt))
        if (self.analysisControl['doInit']):
            self.initStep()
        if (self.analysisControl['solveVstar']):
            self.solveVstar(dt, self.analysisControl['addTransient'])
        if (self.analysisControl['solveP']):
            self.solveP(dt)
        if (self.analysisControl['solveVtilde']):
            self.solveVtilde(dt)
        if (self.analysisControl['solveVenhanced']):
            self.solveVenhanced(dt)
        if (self.analysisControl['updatePosition']):
            self.updateParticleMotion()
        if (self.analysisControl['updateStress']):
            self.updateParticleStress()
    
    def initStep(self):
        # reset nodal mass, momentum, and force
        for nodeList in self.nodes:
            for node in nodeList:
                node.wipe()
            
        # map mass and momentum to nodes
        for cell in self.cells:
            # cell.mapMassToNodes()  # for particle formulation only
            cell.mapMomentumToNodes()
    
    def solveVstar(self, dt, addTransient=False):
        # compute nodal forces from shear
        for i in range(self.nCellsX+1):
            for j in range(self.nCellsY+1):
                self.nodes[i][j].setForce(zeros(2))
                
        for cell in self.cells:
            cell.computeForces(addTransient)
        
        # solve for nodal acceleration a*
        # and update nodal velocity to v*
        for i in range(self.nCellsX+1):
            for j in range(self.nCellsY+1):
                self.nodes[i][j].updateVstar(dt)
    
    def solveP(self, dt):
        ndof = (self.nCellsX+1)*(self.nCellsY+1)
        
        # assemble matrix and force 
        self.FP = zeros(ndof)
        self.KP = zeros((ndof,ndof))
        
        for cell in self.cells:
            ke = cell.GetStiffness()
            fe = cell.GetPforce(dt)
            nodeIndices = cell.getGridCoordinates()
            dof = [ x[0] + x[1]*(self.nCellsX+1)   for x in nodeIndices ]
            
            for i in range(4):
                self.FP[dof[i]] += fe[i]
                for j in range(4):
                    self.KP[dof[i]][dof[j]] += ke[i][j]
                
        # apply boundary conditions
        i = self.nCellsX // 2
        dof = i + self.nCellsY*(self.nCellsX+1)
        self.KP[dof][dof] = 1.0e20
        self.FP[dof] = 0.0
            
        # solve for nodal p
        pressure = solve(self.KP, self.FP)
        
        # assign pressure to nodes
        for i in range(self.nCellsX+1):
            for j in range(self.nCellsY+1):
                dof = i + j*(self.nCellsX+1)
                self.nodes[i][j].setPressure(pressure[dof])
                
        #print(pressure)
        
    def solveVtilde(self, dt):
        for i in range(self.nCellsX+1):
            for j in range(self.nCellsY+1):
                # compute nodal pressure gradient
                if ( i==0 or i==self.nCellsX ):
                    dpx = 0.0
                else:
                    dpx = 0.5*(self.nodes[i+1][j].getPressure()-self.nodes[i-1][j].getPressure())/self.hx
                
                if ( j==0 or j==self.nCellsY ):
                    dpy = 0.0
                else:
                    dpy = 0.5*(self.nodes[i][j+1].getPressure()-self.nodes[i][j-1].getPressure())/self.hy
        
                # update nodal velocity
                dv = -dt/self.rho * array([dpx,dpy])
                self.nodes[i][j].addVelocity(dv)
                
    def solveVenhanced(self, dt):
        for cell in self.cells:
            # initialize the divergence terms in the cell
            cell.SetVelocity()
    
    def updateParticleStress(self):
        pass
    
    def updateParticleMotion(self):
        pass
    
    def findCell(self, x):
        pass
    
    def createParticles(self, n, m):
        for cell in self.cells:
            h = cell.getSize()
            mp = self.rho,h[0]*h[1]/n/m
            
            for i in range(n):
                s = -1. + (2*i+1)/n
                for j in range(j):
                    t = -1. + (2*j+1)/m
                    xl = array([s,t])
                    xp = cell.getGlobal(xl)
                    newParticle = Particle(mp,xp)
                    self.particles.append(newParticle)
                    cell.addParticle(newParticle)
    
    def getTimeStep(self, CFL):
        dt = 1.0e10
        
        for nodeList in self.nodes:
            for node in nodeList:
                vel = node.getVelocity()
                if (abs(vel[0]) > 1.0e-5):
                    dtx = self.hx / abs(vel[0])
                    if (dtx<dt):
                        dt = dtx
                if (abs(vel[1]) > 1.0e-5):
                    dty = self.hy / abs(vel[1])
                    if (dty<dt):
                        dt = dty

        return dt
                