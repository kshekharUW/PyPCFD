'''
Created on Nov 21, 2015

@author: pmackenz
'''

from Domain import *
from ButcherTableau import *
import subprocess

def Main():
    # defne the Reynolds number
    Re = 1000
    
    # set sliding velocity
    velocity = 1.0
    
    # mass density of the fluid
    density = 1000.
    
    # set side-length of the analysis domain
    edgeDomain      = 1.
    # set the number of cells per edge
    numCellsPerEdge = 8
    
    # viscosity of the fluid
    viscosity = density * velocity * edgeDomain / Re
    
    # create an analysis domain
    domain = Domain(edgeDomain,edgeDomain,numCellsPerEdge,numCellsPerEdge)
    domain.createParticles(2, 2)
    domain.setTimeIntegrator(ExplicitEuler())
    
    # configure the analysis type
    doInit         = False
    solveVstar     = True
    solveP         = True
    solveVtilde    = True
    solveVenhanced = False
    updatePosition = True
    updateStress   = False
    addTransient   = True
    
    domain.setAnalysis(doInit,
                       solveVstar,
                       solveP,
                       solveVtilde,
                       solveVenhanced,
                       updatePosition,
                       updateStress,
                       addTransient)
    domain.setParameters(Re, density, velocity)
    domain.setInitialState()
    
    CFL = 1.0
    dt = domain.getTimeStep(CFL)
    
    print(u"CFL=1 equals to \u0394t={:f}".format(dt))


    print(domain)

    # define load history and print interval

    dt1 = 0.5
    target1 = 10.0
    
    dt1 = 0.00125
    dt1 = 0.10000
    target1 = 2.

    dt2 = 0.5
    target2 = 1.0

# ************* don't mess with stuff below *************

    # defining plot settings
    domain.setPlotInterval(dt1)

    # deining output settings
    domain.setWriteInterval(-1)

    # initializing starting time
    time = 0.0
    
    # run first segment
    dt = dt1
    while (time+dt <= target1+0.1*dt):
        time += dt
        domain.runAnalysis(time)
    
    # generate the animation
    subprocess.run('./makeAnim.sh')
    

if __name__ == '__main__':
    Main()
