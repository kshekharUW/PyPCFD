# ====== settings ================

PLOT_MOTIONS = False
PLOT_SINGLE_STEP_TESTS = True
PLOT_MULTI_STEP_TESTS  = True
COLLATE_PLOTS = False

MOTION1 = True
MOTION2 = True
MOTION3 = False
MOTION4 = False

ALGORITHM_EXPLICIT    = True
ALGORITHM_MIDPOINT    = True
ALGORITHM_RUNGE_KUTTA = True

OUTPUT_FILE_TYPE = 'png'

NUM_CELLS = 1

# ====== the test function =======
from MotionPlot import *
from ErrorPlotter import *


def testSingleStep(theMotion):
    if COLLATE_PLOTS:
        errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS)
        if ALGORITHM_EXPLICIT:
            errorPlot.addTestData(LocalConvergenceTest(theMotion, ExplicitEuler(), OUTPUT_FILE_TYPE, NUM_CELLS))

        if ALGORITHM_MIDPOINT:
            errorPlot.addTestData(LocalConvergenceTest(theMotion, MidPointRule(), OUTPUT_FILE_TYPE, NUM_CELLS))

        if ALGORITHM_RUNGE_KUTTA:
            errorPlot.addTestData(LocalConvergenceTest(theMotion, RungeKutta4(), OUTPUT_FILE_TYPE, NUM_CELLS))
        errorPlot.savePlot(theMotion)

    else:
        if ALGORITHM_EXPLICIT:
            numAlg = ExplicitEuler()
            errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS, numAlg, OUTPUT_FILE_TYPE)
            errorPlot.addTestData(
                LocalConvergenceTest(theMotion, numAlg, OUTPUT_FILE_TYPE, nCells=NUM_CELLS))
            errorPlot.savePlot(theMotion)

        if ALGORITHM_MIDPOINT:
            numAlg = MidPointRule()
            errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS, numAlg, OUTPUT_FILE_TYPE)
            errorPlot.addTestData(
                LocalConvergenceTest(theMotion, numAlg, OUTPUT_FILE_TYPE, nCells=NUM_CELLS))
            errorPlot.savePlot(theMotion)

        if ALGORITHM_RUNGE_KUTTA:
            numAlg = RungeKutta4()
            errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS, numAlg, OUTPUT_FILE_TYPE)
            errorPlot.addTestData(
                LocalConvergenceTest(theMotion, numAlg, OUTPUT_FILE_TYPE, nCells=NUM_CELLS))
            errorPlot.savePlot(theMotion)


def testMultipleSteps(theMotion):
    if COLLATE_PLOTS:
        errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS)
        if ALGORITHM_EXPLICIT:
            errorPlot.addTestData(GlobalConvergenceTest(theMotion, ExplicitEuler(), OUTPUT_FILE_TYPE, NUM_CELLS))

        if ALGORITHM_MIDPOINT:
            errorPlot.addTestData(GlobalConvergenceTest(theMotion, MidPointRule(), OUTPUT_FILE_TYPE, NUM_CELLS))

        if ALGORITHM_RUNGE_KUTTA:
            errorPlot.addTestData(GlobalConvergenceTest(theMotion, RungeKutta4(), OUTPUT_FILE_TYPE, NUM_CELLS))
        errorPlot.savePlot(theMotion)

    else:
        if ALGORITHM_EXPLICIT:
            numAlg = ExplicitEuler()
            errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS, numAlg, OUTPUT_FILE_TYPE)
            errorPlot.addTestData(
                GlobalConvergenceTest(theMotion, numAlg, OUTPUT_FILE_TYPE, nCells=NUM_CELLS))
            errorPlot.savePlot(theMotion)

        if ALGORITHM_MIDPOINT:
            numAlg = MidPointRule()
            errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS, numAlg, OUTPUT_FILE_TYPE)
            errorPlot.addTestData(
                GlobalConvergenceTest(theMotion, numAlg, OUTPUT_FILE_TYPE, nCells=NUM_CELLS))
            errorPlot.savePlot(theMotion)

        if ALGORITHM_RUNGE_KUTTA:
            numAlg = RungeKutta4()
            errorPlot = ErrorPlotter(NUM_CELLS, COLLATE_PLOTS, numAlg, OUTPUT_FILE_TYPE)
            errorPlot.addTestData(
                GlobalConvergenceTest(theMotion, numAlg, OUTPUT_FILE_TYPE, nCells=NUM_CELLS))
            errorPlot.savePlot(theMotion)


def plotMotionTraces():

    if MOTION1:
        m = MotionPlot(Motion1())
        m.setMaxTime(4.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m1.png")

    if MOTION2:
        m = MotionPlot(Motion2())
        m.setMaxTime(10.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m2.png")

    if MOTION3:
        m = MotionPlot(Motion3())
        m.setMaxTime(8.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m3.png")

    if MOTION4:
        m = MotionPlot(Motion4())
        m.setMaxTime(20.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m4.png")

        m = MotionPlot(Motion4())
        m.setMaxTime(20.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.1, .5], [0.3, .5], [0.5, .5], [0.7, .5], [0.9, .5]))
        m.exportImage("m4b.png")

def Main():
    if PLOT_SINGLE_STEP_TESTS:

        if MOTION1:
            testSingleStep(Motion1())

        if MOTION2:
            testSingleStep(Motion2())

        if MOTION3:
            testSingleStep(Motion3())

        if MOTION4:
            testSingleStep(Motion4())

    if PLOT_MULTI_STEP_TESTS:

        if MOTION1:
            testMultipleSteps(Motion1())

        if MOTION2:
            testMultipleSteps(Motion2())

        if MOTION3:
            testMultipleSteps(Motion3())

        if MOTION4:
            testMultipleSteps(Motion4())

    if PLOT_MOTIONS:
        plotMotionTraces()


if __name__ == '__main__':
    Main()
