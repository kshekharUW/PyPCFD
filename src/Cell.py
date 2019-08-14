'''
Created on Nov 21, 2015

@author: pmackenz
'''
from numpy import array, dot, cross, outer, tensordot, zeros, ones, sqrt, stack, mat
from _operator import index

class Cell(object):
    '''
    variables:
        self.id     = id
        self.nodes  = []
        self.useEnhanced = False
        self.ux = zeros(4)    # velocity field
        self.uy = zeros(4)    # velocity field
        self.divVa = 0.0
        self.divVb = 0.0
        self.divVc = 0.0
        self.rho  = density
        self.mu   = viscosity
        self.xm   = zeros(2)
        self.size = ones(2)
        self.uHat = array([0.0,0.0])  # enhanced field parameters
        self.fHat = array([0.0,0.0])  # enhanced field forces
        self.mHat = array([0.0,0.0])  # enhanced field mass
        self.setShape(array([0.0,0.0]))
        self.myParticles = []
    
    methods:
        def __init__(self, id)
        def setParameters(self, density, viscosity)
        def setEnhanced(self, useEnhanced=True)
        def addParticle(self, particle)
        def releaseParticles(self)
        def getLocal(self, x)
        def getGlobal(self, xl)
        def setShape(self,xl)
        def SetNodes(self, nds)
        def SetVelocity(self, u)
        def GetVelocity(self, x)
        def GetApparentAccel(self, x)
        def SetPressure(self, p)
        def GetPressure(self, x)
        def GetGradientP(self, x)
        def GetGradientV(self, x)
        def GetStrainRate(self, xl)
        def GetGradientA(self, x)              # returns gradient of acceleration field
        def GetEnhancedStrainRate(self, xl)
        def computeForces(self)                # compute nodal forces from viscous stress and add them to the nodes
        def GetStiffness(self)                 # "stiffness matrix" for pressure calculation
        def GetPforce(self)                    # driving force for pressure
        def contains(self, x)                  # True of global position x is within mapped domain of this cell
        def getSize(self)                      # return (hx, hy)
        def getGridCoordinates(self)
        def mapMassToNodes(self)
        def mapMomentumToNodes(self)
        def GetAcceleration(self, x)
        def getID(self)
    '''

    def __init__(self, id, hx, hy):
        '''
        Constructor
        '''
        self.id     = id
        self.nodes  = []
        self.ux = zeros(4)    # velocity field
        self.uy = zeros(4)    # velocity field
        
        self.ax = zeros(4)    # apparent acceleration field
        self.ay = zeros(4)    # apparent acceleration  field
        
        self.useEnhanced = False
        
        self.divVa = 0.0
        self.divVb = 0.0
        self.divVc = 0.0
        
        self.xm   = zeros(2)
        self.size = array([hx,hy])
        
        self.uHat = array([0.0,0.0])  # enhanced field parameters
        self.fHat = array([0.0,0.0])  # enhanced field forces
        self.mHat = array([0.0,0.0])  # enhanced field mass

        
        self.setShape(array([0.0,0.0]))
        
        self.myParticles = []

        
    
    def __str__(self):
        s = "   cell({}): ({}/{}),({}/{}),({}/{}),({}/{})".format(self.id,
                                                                    *self.nodes[0].getGridCoordinates(),
                                                                    *self.nodes[1].getGridCoordinates(),
                                                                    *self.nodes[2].getGridCoordinates(),
                                                                    *self.nodes[3].getGridCoordinates()
                                                                    )
        return s
    
    def setParameters(self, density, viscosity):
        self.rho = density
        self.mu  = viscosity
    
    def setEnhanced(self, useEnhanced=True):
        self.useEnhanced = useEnhanced
        
    def addParticle(self, particle):
        self.myParticles.append(particle)
        
    def releaseParticles(self):
        listOfReleasedParticles = []
        listOfLocalParticles = []
        
        for p in self.myParticles:
            if ( self.contains(p.position()) ):
                listOfLocalParticles.append(p)
            else:
                listOfReleasedParticles.append(p)
        
        self.myParticles = listOfLocalParticles
        
        return listOfReleasedParticles
    
    def getLocal(self, x):
        xl = 2*(x - self.xm) / self.size
        return xl
    
    def getGlobal(self, xl):
        x = 0.5*xl*self.size + self.xm
        return x
    
    def setShape(self,xl):
        sp = 0.5*(1. + xl[0])
        sm = 0.5*(1. - xl[0])
        tp = 0.5*(1. + xl[1])
        tm = 0.5*(1. - xl[1])
        self.shape   = array([ sm*tm, sp*tm, sp*tp, sm*tp ])
        self.DshapeX = array([ -tm,  tm,  tp, -tp ]) / self.size[0]
        self.DshapeY = array([ -sm, -sp,  sp,  sm ]) / self.size[1]
    
    def SetNodes(self, nds):
        self.nodes = nds
        self.xm = zeros(2)
        for node in nds:
            self.xm += 0.25*node.getPosition()
        
    def GetNodeIndexes(self):
        indexes = []
        for node in self.nodes:
            indexes.append(node.getGridCoordinates())
        return indexes
        
    def SetVelocity(self):
        
        self.ux = zeros(4)
        self.uy = zeros(4)
        
        for i in range(4):
            vel = self.nodes[i].getVelocity()
            self.ux[i] = vel[0]
            self.uy[i] = vel[1]
        
        self.divVa =  0.5*(-self.ux[0] + self.ux[1] + self.ux[2] - self.ux[3]) / self.size[0]
        self.divVa += 0.5*(-self.uy[0] - self.uy[1] + self.uy[2] + self.uy[3]) / self.size[1]
        self.divVb =  0.5*(self.uy[0]-self.uy[1]+self.uy[2]-self.uy[3]) / self.size[1]
        self.divVc =  0.5*(self.ux[0]-self.ux[1]+self.ux[2]-self.ux[3]) / self.size[0]

    # def SetVelocity(self, ux, uy):
    #      for i in range(4):
    #         self.nodes[i].SetVelocity(array([ux]))

        
    def GetVelocity(self, x):
        xl = self.getLocal(x)
        self.setShape(xl)
        vel = array([dot(self.shape, self.ux), dot(self.shape, self.uy)])
        
        if (self.useEnhanced):
            # add the enhanced velocity field
            dvx = 0.5 * self.divVb * (1. - xl[0]*xl[0]) 
            dvy = 0.5 * self.divVc * (1. - xl[1]*xl[1])
            vel += array([dvx, dvy])
            
        return vel

    def GetApparentAccel(self, x):
        xl = self.getLocal(x)
        self.setShape(xl)
        accel = array([dot(self.shape, self.ax), dot(self.shape, self.ay)])
            
        return accel

    def GetAcceleration(self, x):
        xl = self.getLocal(x)
        self.setShape(xl)
        ax = zeros((4,1))
        ay = zeros((4,1))
        for i in range(4):
            nodalforce = self.nodes[i].getForce()
            ax[i] = nodalforce[0]/self.nodes[i].getMass()
            ay[i] = nodalforce[1]/self.nodes[i].getMass()
        accn = array([dot(self.shape, ax), dot(self.shape, ay)])
            
        return accn

    def SetPressure(self, p):
        self.p = p

    def GetPressure(self, x):
        xl = self.getLocal(x)
        self.setShape(xl)
        return dot(self.shape, self.p)
    
    def GetGradientP(self, x):
        xl = self.getLocal(x)
        self.setShape(xl)
        
        return array([dot(self.DshapeX,self.p), dot(self.DshapeY,self.p) ])

    def GetGradientV(self, x):
        xl = self.getLocal(x)
        self.setShape(xl)

        dxu = dot(self.DshapeX, self.ux)
        dyu = dot(self.DshapeY, self.ux)
        dxv = dot(self.DshapeX, self.uy)
        dyv = dot(self.DshapeY, self.uy)

        return array([[dxu, dyu],[dxv, dyv]])

    # HACK to test deformation Gradient
    def GetGradientA(self, x): 
        xl = self.getLocal(x)
        self.setShape(xl)
        
        self.ax = zeros(4)
        self.ay = zeros(4)
        
        for i in range(4):
            accel = self.nodes[i].getApparentAccel()
            self.ax[i] = accel[0]
            self.ay[i] = accel[1]
        
        dxax = dot(self.DshapeX, self.ax)
        dyax = dot(self.DshapeY, self.ax)
        dxay = dot(self.DshapeX, self.ay)
        dyay = dot(self.DshapeY, self.ay)
            
        return array([[dxax, dyax],[dxay, dyay]])
    
    def GetStrainRate(self, xl):
        self.setShape(xl)
        
        dxu = dot(self.DshapeX, self.ux)
        dyu = dot(self.DshapeY, self.ux)
        dxv = dot(self.DshapeX, self.uy)
        dyv = dot(self.DshapeY, self.uy)

        #return array([dxu, dyv, dyu+dxv])
        
        dd = (dxu + dyv) / 3.
        return array([dxu-dd, dyv-dd, dyu+dxv])
    
    def GetEnhancedStrainRate(self, xl):
        
        s = xl[0]
        t = xl[1]
        
        # this is the full rate of deformation tensor
        d = [ -2.*self.divVb*s/self.size[0], 
              -2.*self.divVc*t/self.size[1], 
              0.0 ]

        ## this is the deviatoric portion
        #d = [ (self.divVc*t/self.size[1] - 2.*self.divVb*s/self.size[0])/1.5, 
        #      (self.divVb*s/self.size[0] - 2.*self.divVc*t/self.size[1])/1.5, 
        #      0.0 ]

        return array(d)
    
    def computeForces(self, addTransient=False):
        gpts = [ -1./sqrt(3.), 1./sqrt(3.) ]
        w = self.size[0]*self.size[1]/4.
        
        self.SetVelocity()   # this initializes nodal velocities
        
        forces = zeros([4,2])
        
        for s in gpts:
            for t in gpts:
                xl = array([s,t])
                dh   = self.GetStrainRate(xl)
                
                if (self.useEnhanced):
                    denh = self.GetEnhancedStrainRate(xl)
                else:
                    denh = zeros(3)
                
                d11 = w* 2.0*self.mu * ( dh[0] + denh[0] )
                d22 = w* 2.0*self.mu * ( dh[1] + denh[1] )
                d12 = w*     self.mu * ( dh[2] + denh[2] )
                d21 = d12
                
                dfx = d11*self.DshapeX + d12*self.DshapeY
                dfy = d21*self.DshapeX + d22*self.DshapeY
                
                forces -= stack((dfx,dfy),-1)
                
                if (addTransient):
                    
                    aTransient = zeros(2)
        
                    dxu = dot(self.DshapeX, self.ux)
                    dyu = dot(self.DshapeY, self.ux)
                    dxv = dot(self.DshapeX, self.uy)
                    dyv = dot(self.DshapeY, self.uy)
                    
                    # add  w . (grad v) . v
                    vx = dot(self.shape, self.ux)
                    vy = dot(self.shape, self.uy)
                    # standard tensor (single) dot product
                    aTransient[0] = dxu * vx + dyu * vy 
                    aTransient[1] = dxv * vx + dyv * vy 
                    
                    fTransient = w * self.rho * tensordot(self.shape, aTransient, axes=0)  # tensor product
                    
                    forces -= fTransient
                
                
        for i in range(4):
            self.nodes[i].addForce(forces[i])
            
        
    def GetStiffness(self):
        gpts = [ -1./sqrt(3.), 1./sqrt(3.) ]
        w = self.size[0]*self.size[1]/4.
        
        Ke = zeros((4,4))
        
        for s in gpts:
            for t in gpts:
                xl = array([s,t])
                self.setShape(xl)
                
                B = stack((self.DshapeX,self.DshapeY))
                
                Ke += w*tensordot(B, B, ([0,0]))
        
        return Ke
    
    def GetPforce(self,dt):
        gpts = [ -1./sqrt(3.), 1./sqrt(3.) ]
        w = self.rho*self.size[0]*self.size[1]/4./dt
        
        self.SetVelocity()
        
        Fe = zeros(4)
        
        for s in gpts:
            for t in gpts:
                xl = array([s,t])
                self.setShape(xl)
                
                divV = self.divVa + self.divVb*xl[0] + self.divVc*xl[1]
                
                Fe += -w*self.shape*divV
        
        return Fe
    
    def contains(self, x):
        xl = self.getLocal(x)
        state = True
        if (xl[0]<-1. or xl[0]>+1. or xl[1]<-1. or xl[1]>+1.):
            state = False
            
        return state
    
    def getSize(self):
        return self.size
    
    def getGridCoordinates(self):
        coords = []
        for node in self.nodes:
            coords.append(node.getGridCoordinates())
        return coords
    
    def mapMassToNodes(self):
        gpts = [ -1./sqrt(3.), 1./sqrt(3.) ]
        w = self.rho * self.size[0]*self.size[1]/4.
        
        mass     = zeros(4)
        
        for s in gpts:
            for t in gpts:
                xl = array([s,t])
                self.setShape(xl)
                mass += w*self.shape
                
        for i in range(4):
            self.nodes[i].addMass(mass[i])
            
    def mapMomentumToNodes(self):
        nodalV = zeros((2,4))
        for i in range(4):
            nodalV[:,i] = self.nodes[i].getVelocity()
            
        gpts = [ -1./sqrt(3.), 1./sqrt(3.) ]
        w = self.rho * self.size[0]*self.size[1]/4.
        
        momentum = zeros((2,4))
        
        for s in gpts:
            for t in gpts:
                xl = array([s,t])
                self.setShape(xl)
                vel = tensordot(nodalV, self.shape, ([1,0])) 
                momentum += w*self.shape*vel
                
        for i in range(4):
            self.nodes[i].addMomentum(momentum[i])

    def getID(self):
        return self.id



        
        
            
        