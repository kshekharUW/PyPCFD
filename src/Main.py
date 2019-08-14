'''
Created on Nov 21, 2015

@author: pmackenz
'''

from Domain import *
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
    
    # configure the analysis type
    doInit         = False
    solveVstar     = True
    solveP         = True
    solveVtilde    = True
    solveVenhanced = False
    updatePosition = False
    updateStress   = False
    addTransient   = False
    plotFigures    = True
    writeOutput    = False
    
    domain.setAnalysis(doInit,
                       solveVstar,
		       solveP,
		       solveVtilde,
		       solveVenhanced,
		       updatePosition,
		       updateStress,
		       addTransient,
		       plotFigures,
		       writeOutput)
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
    target1 = 10.0

    dt2 = 0.5
    target2 = 1.0

# ************* don't mess with stuff below *************

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
