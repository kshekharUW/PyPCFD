from Domain import *
import matplotlib.pyplot as plt
from math import log


def Main():
    # configure the analysis type
    doInit         = False
    solveVstar     = False
    solveP         = False
    solveVtilde    = False
    solveVenhanced = False
    updatePosition = True
    updateStress   = False
    addTransient   = False
    plotFigures    = True
    writeOutput    = False

    errorList = []
    dtList = []
    dt = 1.0
    while ( dt>1.0e-8 ):
        dtList.append(dt)
        domain = Domain(width=2, height=2, nCellsX=1, nCellsY=1)
        domain.setAnalysis(doInit, solveVstar, solveP, solveVtilde, solveVenhanced, updatePosition, updateStress, addTransient, plotFigures, writeOutput)
        #domain.setState(0)
        # domain.plotData()
        domain.setState(dt)
        errorList.append(domain.updateParticleMotion(dt)[0]) # returns error for all particles. Use [0] for first [1] for second etc..
        print('dt = {:.2E}, error = {:.3E}'.format(dt, errorList[-1]))
        if (errorList[-1]<1.e-16):
            break
        dt /= 10.0
        # domain.plotData()

    plt.figure(1)
    plt.loglog(dtList, errorList, 'r-o')
    #plt.scatter(dtList, errorList)  # not necessary
    
    x = array([dtList[0],dtList[-1]])
    y = array([errorList[0], errorList[0]*(dtList[-1] / dtList[0])**(2.)])
    plt.loglog(x,y,'b--')
    
    x = array([dtList[0],dtList[-1]])
    y = array([errorList[0], errorList[0]*(dtList[-1] / dtList[0])**(3.)])
    plt.loglog(x,y,'g--')
    
    x = array([dtList[0],dtList[-1]])
    y = array([errorList[0], errorList[0]*(dtList[-1] / dtList[0])**(4.)])
    plt.loglog(x,y,'m--')
    
    x = array([dtList[0],dtList[-1]])
    y = array([errorList[0], errorList[0]*(dtList[-1] / dtList[0])**(5.)])
    plt.loglog(x,y,'y--')
    
    plt.xlabel('$\Delta t$ (s)')
    plt.ylabel('$|| F_{numerical} - F_{analytical} ||_{2}$')
    
    plt.legend(("simulation", "2nd order", "3rd order", "4th order", "5th order"))
    
    plt.grid(True)
    plt.axis('tight')
    plt.savefig("convergence.pdf", bbox_inches='tight')
    plt.savefig("convergence.png", bbox_inches='tight')

    slope = log(errorList[0] / errorList[-1]) / log(dtList[0] / dtList[-1])
    print('slope = {:.2E}'.format(slope))


if __name__ == '__main__':
    Main()