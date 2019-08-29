import matplotlib.pyplot as plt
from math import log
from LocalConvergenceTest import *
from GlobalConvergenceTest import *

class ErrorPlotter(object):

    def __init__(self):

        self.xLabel = "'$\Delta t$ (s)'"
        self.folderName = "Single_Step"

        self.figF, self.axF = plt.subplots()
        self.figP, self.axP = plt.subplots()

        self.createFolders()

    def createFolders(self):
        if not os.path.isdir("images"):
            os.mkdir("images")

        if not os.path.isdir(os.path.join("images", "Single_Step")):
            os.mkdir(os.path.join("images", "Single_Step"))

        if not os.path.isdir(os.path.join("images", "Multi_Step")):
            os.mkdir(os.path.join("images", "Multi_Step"))

    def getNumAlgLineStyle(self, numAlg):
        if isinstance(numAlg, ExplicitEuler):
            return "k-o"
        elif isinstance(numAlg, MidPointRule):
            return "r-s"
        elif isinstance(numAlg, RungeKutta4):
            return "b-^"

    def addTestData(self, test):
        if isinstance(test, GlobalConvergenceTest):
            self.xLabel = '$N$'
            self.folderName = "Multi_Step"

        data = test.getErrors()
        self.plotPositionErrors(data[0], data[1], test.getNumAlg(), test.getMotion())
        self.plotFErrors(data[0], data[2], test.getNumAlg(), test.getMotion())

    def plotPositionErrors(self, xList, positionErrors, numAlg, motion):
        self.axP.loglog(xList, positionErrors, self.getNumAlgLineStyle(numAlg), linewidth=2, label=numAlg)

        slope = log(positionErrors[0] / positionErrors[-1]) / log(xList[0] / xList[-1])
        print('{} {} Position convergence slope = {:.2f}'.format(numAlg, motion, slope))

    def plotFErrors(self, xList, FErrors, numAlg, motion):
        self.axF.loglog(xList, FErrors, self.getNumAlgLineStyle(numAlg), linewidth=2, label=numAlg)

        slope = log(FErrors[0] / FErrors[-1]) / log(xList[0] / xList[-1])
        print('{} {} F convergence slope = {:.2f}'.format(numAlg, motion, slope))

    def savePlot(self, motion):
        self.axP.set_xlabel(self.xLabel)
        self.axP.set_ylabel('$|| x_{numerical} - x_{analytical} ||_{2}$')
        self.axP.legend(loc="best")
        self.axP.grid(True)
        fileName = "{}_Position_convergence.{}".format(motion, 'png')
        fileNameWithPath = os.path.join("images", self.folderName, fileName)
        self.figP.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')


        self.axF.set_xlabel(self.xLabel)
        self.axF.set_ylabel('$|| F_{numerical} - F_{analytical} ||_{2}$')
        self.axF.legend(loc="best")
        self.axF.grid(True)
        fileName = "{}_F_convergence.{}".format(motion, 'png')
        fileNameWithPath = os.path.join("images", self.folderName, fileName)
        self.figF.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')

