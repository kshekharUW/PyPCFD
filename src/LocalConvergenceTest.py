from math import floor
from numpy.linalg import norm

from Domain import *
from Motion import *
from ButcherTableau import *


class LocalConvergenceTest(object):

    def __init__(self, motion, algorithm, fileType='png', nCells=1):
        self.numAlgorithm = algorithm
        self.motion = motion
        self.fileType = fileType
        if nCells >= 1.:
            self.nCells= floor(nCells)
        else:
            self.nCells = 1

        # Error storage containers
        self.Ferrors = []
        self.positionErrors = []
        self.dtList = []

        self.runCase(self.numAlgorithm, self.motion)

    def getErrors(self):
        return self.dtList, self.positionErrors, self.Ferrors

    def getNumAlg(self):
        return self.numAlgorithm

    def getMotion(self):
        return self.motion

    def runCase(self, numAlg, motion):
        # configure the analysis type
        self.doInit = False
        self.solveVstar = False
        self.solveP = False
        self.solveVtilde = False
        self.solveVenhanced = False
        self.updatePosition = True
        self.updateStress = False
        self.addTransient = False

        dt = 1.0
        while (dt > 1.0e-10):
            self.dtList.append(dt)
            width = 1.
            height = 1.
            domain = Domain(width=height, height=height, nCellsX=self.nCells, nCellsY=self.nCells)
            domain.setMotion(motion)
            domain.setTimeIntegrator(numAlg)

            domain.setAnalysis(self.doInit, self.solveVstar, self.solveP,
                               self.solveVtilde, self.solveVenhanced,
                               self.updatePosition, self.updateStress,
                               self.addTransient)
            domain.createParticleAtX(1.0, array([width/2.,height/10.]))

            domain.setPlotInterval(dt)        # plot at the end of each time step
            domain.setWriteInterval(-1)       # no recorder output

            # you need to set the velocity field to the initial velocity field
            # (or to any fixed time throughout the test !!!)
            domain.setState(0)
            x0 = domain.getParticles()[0].position() # save original position of particle

            # update particle
            domain.updateParticleMotion(dt)
            # calculate and store errors from updated particle position
            particle = domain.getParticles()[0]
            posError = norm( particle.position() - motion.getAnalyticalPosition(x0, dt) )
            FError = norm(particle.getDeformationGradient() - motion.getAnalyticalF(x0, dt) )
            self.Ferrors.append(FError)
            self.positionErrors.append(posError)

            mask = 'dt = {:.2E}, Position error = {:.3E}, F error = {:.3E}'
            print(mask.format(dt, self.positionErrors[-1], self.Ferrors[-1]))
            if self.Ferrors[-1] < 1.e-15 or self.positionErrors[-1] < 1.e-15:
                break
            # dt /= 10.
            dt /= 2.  # allows for a more precise identification of numeric truncation error
