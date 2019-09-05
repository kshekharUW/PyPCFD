from math import floor
from numpy.linalg import norm

from Domain import *
from Motion import *
from ButcherTableau import *


class GlobalConvergenceTest(object):

    def __init__(self, motion, algorithm, fileType='png', nCells=1):
        self.numAlgorithm = algorithm
        self.motion = motion
        self.fileType = fileType
        if nCells >= 1.:
            self.nCells= floor(nCells)
        else:
            self.nCells = 1

        # configure the analysis type
        self.doInit = False
        self.solveVstar = False
        self.solveP = False
        self.solveVtilde = False
        self.solveVenhanced = False
        self.updatePosition = True
        self.updateStress = False
        self.addTransient = False

        # Domain parameters
        self.domainHeight = 1.0
        self.domainWidth = 1.0

        # Error storage containers
        self.Ferrors = []
        self.positionErrors = []
        self.NList = []

        self.runCase(self.numAlgorithm, self.motion)

    def getErrors(self):
        return self.NList, self.positionErrors, self.Ferrors

    def getNumAlg(self):
        return self.numAlgorithm

    def getMotion(self):
        return self.motion

    def runCase(self, numAlg, motion):
        maxTime = 1.0
        N = 1
        dt = maxTime/N
        while (N<=1000):
            self.NList.append(N)

            width = 1.
            height = 1.
            domain = Domain(width=width, height=height, nCellsX=self.nCells, nCellsY=self.nCells)
            domain.setMotion(motion)
            domain.setTimeIntegrator(numAlg)

            domain.setAnalysis(self.doInit, self.solveVstar, self.solveP,
                               self.solveVtilde, self.solveVenhanced,
                               self.updatePosition, self.updateStress,
                               self.addTransient)
            domain.createParticleAtX(1.0, array([width / 2., height / 10.]))

            domain.setPlotInterval(maxTime)   # plot only at the end
            domain.setWriteInterval(-1)       # no recorder output

            # Set the velocity field to the initial velocity field
            x0 = domain.getParticles()[0].position()  # save original position of particle for comparison later

            # update particle
            for j in range(N):
                domain.setState(j*dt)
                domain.updateParticleMotion(dt)

            # calculate and store errors from updated particle position
            particle = domain.getParticles()[0]
            posError = norm( particle.position() - motion.getAnalyticalPosition(x0, dt*N))
            FError = norm(particle.getDeformationGradient() - motion.getAnalyticalF(x0, dt*N))
            self.Ferrors.append(FError)
            self.positionErrors.append(posError)

            mask = 'N = {}, dt = {:.2E}, Position error = {:.3E}, F error = {:.3E}'
            print(mask.format(N, dt, self.positionErrors[-1], self.Ferrors[-1]))
            if (self.Ferrors[-1] < 1.e-14 or self.positionErrors[-1] < 1.e-14):
                break
            N *= 10
            dt = maxTime/N
