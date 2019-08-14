'''
Created on Nov 21, 2015

@author: pmackenz
'''
from numpy import array, zeros, ones
from numpy.linalg import norm
from sympy.physics.units.dimensions import mass
from scipy.linalg import expm
from math import pi

class Node(object):
    '''
    variables:
        self.id = id
        self.pos = array([X,Y])
        self.force = zeros(2)
        self.mass = 0.0
        self.momentum = zeros(2)
        self.pressure
        self.aStar  = zeros(2)
        self.aTilde = zeros(2)
        self.appAccel = zeros(2)
        self.ahat   = zeros(2)  # should always remain zero for interpolation purposes
        self.lastX = self.pos   # last converged position
        self.lastV = zeros(2)   # last converged velocity
        self.fixety = dict()
        self.gridCoords = (i,j)
    
    methods:
        def __init__(self, id, X,Y)
        def __str__(self)
        def wipe(self)
        def setGridCoordinates(self, i,j)
        def getGridCoordinates(self)
        def getPosition(self)
        def setMomentum(self, p)
        def addMomentum(self, p)
        def getMomentum(self)
        def setMass(self, m)
        def addMass(self, m)
        def getMass(self)
        def setVelocity(self, v)
        def addVelocity(self, dv)
        def getVelocity(self)
        def setApparentAccel(self, a)
        def getApparentAccel(self)
        def setPressure(self, p)
        def getPressure(self)
        def setForce(self, F)
        def addForce(self, F)
        def getForce(self)
        def getFixeties(self)
        def fixDOF(self, i, val=0.0)
        def updateVstar(self)
        def updateV(self, v)
    '''


    def __init__(self, id, X,Y):
        '''
        Constructor
        '''
        self.id = id
        self.gridCoords = ()
        
        self.pos = array([X,Y])
        
        self.force = zeros(2)
        
        self.mass = 0.0
        self.momentum = zeros(2)
        self.appAccel = zeros(2)
        self.pressure = 0.0
        
        self.aStar  = zeros(2)
        self.aTilde = zeros(2)
        self.ahat   = zeros(2)  # should always remain zero for interpolation purposes
        
        self.lastX = self.pos   # last converged position
        self.lastV = zeros(2)   # last converged velocity
        
        self.fixety = dict()
    
    def __str__(self):
        s = "   node({}/{}):  x=[{},{}], mass={}, p=[{},{}], v=[{},{}]".format(*self.gridCoords,
                                                                    *self.pos,
                                                                    self.mass,
                                                                    *self.momentum,
                                                                    *(self.momentum/self.mass))
        return s
    
    def wipe(self):
        self.mass = 0.0
        self.momentum = 0.0
        self.force = zeros(2)
        
    def setGridCoordinates(self, i,j):
        self.gridCoords = (i,j)
        
    def getGridCoordinates(self):
        return self.gridCoords
    
    def getPosition(self):
        return self.pos
                           
    def setMomentum(self, p):
        self.momentum = p
        
    def addMomentum(self, p):
        self.momentum += p
        
    def getMomentum(self):
        return self.momentum
    
    def setMass(self, m):
        self.mass = m
    
    def addMass(self, m):
        self.mass += m

    def getMass(self):
        return self.mass
    
    def setVelocity(self, v):
        self.momentum = self.mass*v
        
    def addVelocity(self, dv):
        # check for boundary conditions !!!!
        if (not 0 in self.fixety):
            self.momentum[0] += self.mass*dv[0]
        if (not 1 in self.fixety):
            self.momentum[1] += self.mass*dv[1]
        
    def getVelocity(self):
        if (self.mass > 0.0):
            return self.momentum/self.mass
        else:
            print("NO mass at node")
            raise
    
    def setApparentAccel(self, a):
        self.appAccel = a
        
    def getApparentAccel(self):
        return self.appAccel
        
    def setPressure(self, p):
        self.pressure = p
        
    def getPressure(self):
        return self.pressure
    
    def setForce(self, F):
        self.force = F
        
    def addForce(self, F):
        self.force += F
        
    def getForce(self):
        return self.force
        
    def getFixeties(self):
        fixeties = []
        for i in self.fixety.keys():
            fixeties.append((i,fixeties[i]))
        return fixeties
    
    def fixDOF(self, dof, val=0.0):
        self.fixety[dof] = val
       
    def updateVstar(self, dt):
        # apply boundary condition
        for dof in self.fixety.keys():
            val = self.fixety[dof]
            self.force[dof] = (self.mass * val - self.momentum[dof])/dt
        
        # update velocity
        self.aStar = self.force / self.mass
        self.addVelocity(self.aStar * dt)
        
    def updateV(self, v):
        pass
        
   

