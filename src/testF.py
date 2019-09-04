# ====== settings ================

PLOT_MOTIONS = False

PLOT_SINGLE_STEP_TESTS = True

PLOT_MULTI_STEP_TESTS  = True

MOTION1 = True
MOTION2 = True
MOTION3 = False
MOTION4 = False

ALGORITHM_EXPLICIT    = True
ALGORITHM_MIDPOINT    = True
ALGORITHM_RUNGE_KUTTA = True

OUTPUT_FILE_TYPE = 'png'

NUM_CELLS = 8

# ====== the test function =======
from MotionPlot import *
from ErrorPlotter import *

def Main():
    fileType = OUTPUT_FILE_TYPE

    if PLOT_SINGLE_STEP_TESTS:

        if MOTION1:
            localErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion1(), ExplicitEuler(), fileType, NUM_CELLS))

            if ALGORITHM_MIDPOINT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion1(), MidPointRule(), fileType, NUM_CELLS))

            if ALGORITHM_RUNGE_KUTTA:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion1(), RungeKutta4(), fileType, NUM_CELLS))
            localErrorPlot.savePlot(Motion1())

        if MOTION2:
            localErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion2(), ExplicitEuler(), fileType, NUM_CELLS))

            if ALGORITHM_MIDPOINT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion2(), MidPointRule(), fileType, NUM_CELLS))

            if ALGORITHM_RUNGE_KUTTA:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion2(), RungeKutta4(), fileType, NUM_CELLS))
            localErrorPlot.savePlot(Motion2())

        if MOTION3:
            localErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion3(), ExplicitEuler(), fileType, NUM_CELLS))

            if ALGORITHM_MIDPOINT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion3(), MidPointRule(), fileType, NUM_CELLS))

            if ALGORITHM_RUNGE_KUTTA:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion3(), RungeKutta4(), fileType, NUM_CELLS))
            localErrorPlot.savePlot(Motion3())

        if MOTION4:
            localErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion4(), ExplicitEuler(), fileType, NUM_CELLS))

            if ALGORITHM_MIDPOINT:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion4(), MidPointRule(), fileType, NUM_CELLS))

            if ALGORITHM_RUNGE_KUTTA:
                localErrorPlot.addTestData(LocalConvergenceTest(Motion4(), RungeKutta4(), fileType, NUM_CELLS))
            localErrorPlot.savePlot(Motion4())


    if PLOT_MULTI_STEP_TESTS:

        if MOTION1:
            globalErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion1(), ExplicitEuler(), fileType))

            if ALGORITHM_MIDPOINT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion1(), MidPointRule(), fileType))

            if ALGORITHM_RUNGE_KUTTA:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion1(), RungeKutta4(), fileType))

            globalErrorPlot.savePlot(Motion1())

        if MOTION2:
            globalErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion2(), ExplicitEuler(), fileType))

            if ALGORITHM_MIDPOINT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion2(), MidPointRule(), fileType))

            if ALGORITHM_RUNGE_KUTTA:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion2(), RungeKutta4(), fileType))

            globalErrorPlot.savePlot(Motion2())

        if MOTION3:
            globalErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion3(), ExplicitEuler(), fileType))

            if ALGORITHM_MIDPOINT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion3(), MidPointRule(), fileType))

            if ALGORITHM_RUNGE_KUTTA:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion3(), RungeKutta4(), fileType))

            globalErrorPlot.savePlot(Motion3())

        if MOTION4:
            globalErrorPlot = ErrorPlotter(NUM_CELLS)
            if ALGORITHM_EXPLICIT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion4(), ExplicitEuler(), fileType))

            if ALGORITHM_MIDPOINT:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion4(), MidPointRule(), fileType))

            if ALGORITHM_RUNGE_KUTTA:
                globalErrorPlot.addTestData(GlobalConvergenceTest(Motion4(), RungeKutta4(), fileType))

            globalErrorPlot.savePlot(Motion4())


    if PLOT_MOTIONS:

        if MOTION1:
            m = MotionPlot(Motion1())
            m.setMaxTime(4.)
            m.setPointsPerSecond(10)
            m.setTracers( ([0.5, .1],[0.5, .3],[0.5, .5],[0.5, .7],[0.5, .9]) )
            m.exportImage("m1.png")

        if MOTION2:
            m = MotionPlot(Motion2())
            m.setMaxTime(10.)
            m.setPointsPerSecond(10)
            m.setTracers( ([0.5, .1],[0.5, .3],[0.5, .5],[0.5, .7],[0.5, .9]) )
            m.exportImage("m2.png")

        if MOTION3:
            m = MotionPlot(Motion3())
            m.setMaxTime(8.)
            m.setPointsPerSecond(10)
            m.setTracers( ([0.5, .1],[0.5, .3],[0.5, .5],[0.5, .7],[0.5, .9]) )
            m.exportImage("m3.png")

        if MOTION4:
            m = MotionPlot(Motion4())
            m.setMaxTime(20.)
            m.setPointsPerSecond(10)
            m.setTracers( ([0.5, .1],[0.5, .3],[0.5, .5],[0.5, .7],[0.5, .9]) )
            m.exportImage("m4.png")

            m = MotionPlot(Motion4())
            m.setMaxTime(20.)
            m.setPointsPerSecond(10)
            m.setTracers( ([0.1, .5],[0.3, .5],[0.5, .5],[0.7, .5],[0.9, .5]) )
            m.exportImage("m4b.png")


if __name__ == '__main__':
    Main()
