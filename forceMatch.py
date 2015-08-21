#!/usr/bin/python2.7

__author__ = 'Greg Poisson'

import sys
import numpy
import time
import math
import forceMatchInit
import forceMatchAtomGroup
import forceMatchProcesses
import forceMatchPlot



class fileVariables:
    # FILE VARIABLES
    configFile = "nacl.config"
    outFileName = "forceMatch_out.dat"
    psf = None
    forceDcd = None
    coordDcd = None
    param = None
    force = None
    coord = None
    temperature = None
    dims = None
    groups = None
    groupsCoord = None
    groupsForce = None
    ionsCoord = None
    ionsForce = None
    atomGroups = None

    debug = False
    grouping = False
    groupingByResName = False
    junkCounter = 0     # counter used for debugging


    # DEFUALT GLOBAL VARIABLES
    rMin = 0
    rMax = 25
    binSize = 0.1
    binCount = 0
    elecPermittivity = 8.854e-12      # electrical permittivity of vacuum C^2/Jm
    boltzmann = 1.9872041e-3          # boltzmann constant in kcal/(K mol)
    exTs = 0            # example time step used for log purposes

    # PLOT DATA VARIABLE
    plots = []
    plotOrder = []

    titles = [
        "Mean Force",
        "Potential of Mean Force",
        "Radial Distribution Frequency",
        "Free Energy",
        "Probability Density Distribution of Distance Frequency",
        "Frequency Distribution of Particle Distance"
    ]

    shortTitles = [
        "MF",
        "PMF",
        "RDF",
        "FE",
        "PROB",
        "FREQ"
    ]

    yLabels = [
        "Avg Force",
        "Free Energy",
        "g(r)",
        "Free Energy \n(kcal/mol)",
        "Probability",
        "Occurrences"
    ]


    '''
        DATA ORGANIZATION:
            Plots -- List of name / dataset pairs
            Data -- Ordered list of data sets. Real-time data populates the front of the list,
                    while post-process data populates the end of the list.
                    Plots may be drawn in any order, but they are built and stored as follows:
            Data set order:

                (Active measurements:)
                [0] -- Mean Force

                (Post-process measurements:)
                [-5] -- Potential of Mean Force
                [-4] -- Radial Distribution Frequency
                [-3] -- Free Energy
                [-2] -- Probability Density Distribution of Distance Frequency
                [-1] -- Frequency of Particle Distance

                EX: Plots-  0       1       2       3       4       5
                         SOD-SOD   Data  SOD-CLA   Data  CLA-CLA   Data
                                 [MF]            [MF]            [TF]
                                 [PMF]           [PMF]           [PMF]
                                 [RDF]           [RDF]           [RDF]
                                 [FE]            [FE]            [FE]
                                 [Prob Dist]     [Prob Dist]     [Prob Dist]
                                 [Freq]          [Freq]          [Freq]

    '''

# Identifies current particle pair, returns True if pair is unique (doesn't yet exist in plots array), False otherwise
def pairIsUnique(a, b):
    global variables
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in variables.plots:
        return False
    elif pairFlipped in variables.plots:
        return False
    else:
        return True

# Builds a blank data set in which to accumulate data
def buildBlankDataSet():

    totForce = numpy.zeros(variables.binCount)
    integF = numpy.zeros(variables.binCount)

    freeEnergy = numpy.zeros(variables.binCount)
    rdf = numpy.zeros(variables.binCount)
    probDensity = numpy.zeros(variables.binCount) # Probability density of data counts
    dataCounts = numpy.zeros(variables.binCount)  # Keep track of the number of measurements made
                                        # in each bin, for averaging purposes

    dataSet = [totForce, integF, freeEnergy, rdf, probDensity, dataCounts]
    return dataSet

# Generates a LAMMPS input file and particle-pair .dat files
def giveOutFiles(outFile):
    top = True
    for p in range(len(variables.plots)):
        if p % 2 == 0:
            outFile.write("{}\nN\t{}\n\n".format(variables.plots[p], variables.binCount))
            outFile2Name = list(variables.plots[p])
            for c in range(len(outFile2Name)):
                if outFile2Name[c] == " ":
                    outFile2Name[c] = "_"
                    top = True
            outFile2Name = "".join(outFile2Name)
            outFile2 = open("{}.dat".format(outFile2Name), 'w')
        else:
            count = 1
            for d in range(variables.binCount):
                if not (math.isnan(variables.plots[p][3][d])):
                    pmf = "{0:.06f}".format(float(variables.plots[p][-5][d]))
                    mf = "{0:.06f}".format(float(variables.plots[p][0][d]))
                    outFile.write("\t{}\t\t{}\t\t{}\t\t{}\n".format(count, (variables.binSize * d), pmf, mf))
                    count += 1
                    for s in range(len(variables.plots[p])):
                        if top:
                            outFile2.write("#\t\t{}\t\t\t\t{}\t\t\t\t{}\t\t\t\t{}\t\t\t\t{}\t\t\t{}\n{}".format(variables.shortTitles[0], variables.shortTitles[1], variables.shortTitles[2], variables.shortTitles[3], variables.shortTitles[4], variables.shortTitles[5], d * variables.binSize))
                            top = False
                        field = "{0:.6f}".format(float(variables.plots[p][s][d]))
                        outFile2.write("\t\t{}".format(field))
                    if ((d+1) * variables.binSize) < variables.rMax:
                        outFile2.write("\n{}".format((d+1) * variables.binSize))
    outFile.close()

cF = None
# if a config file is given as a parameter, assign its name to configFile
if len(sys.argv) > 1:
    cF = sys.argv[1]

# print program banner
print "\n\n\t**** forceMatch.py ****\n"

# initialize variables
variables = fileVariables
if cF is not None:
    variables.configFile = cF

# run initiation scripts
start = time.time()                     # execute a program runtime timer
forceMatchInit.getFiles(variables)      # retrieves file names and grouping settings from config

# determine requested data
forceMatchProcesses.getPlotOrder(variables.configFile, variables)

numTimesteps = len(variables.coord.trajectory)

testFrames = 20

for ts in variables.coord.trajectory:
    tsf = variables.force.trajectory[ts.frame-1]
    if ts.frame < numTimesteps:
        # course grain out solvent data, evaluate grouping
        if ts.frame == 1:
            forceMatchAtomGroup.groupAtoms(variables)
            forceMatchProcesses.executeProcesses(variables, firstFrame=True)
            print "Process MD frames..."

        else:
            forceMatchAtomGroup.updateGroups(variables)
            forceMatchProcesses.executeProcesses(variables)
    sys.stdout.write("Processing MD frames: {0:.2f}% Complete\r".format((float(ts.frame + 1) / float(numTimesteps)) * 100))
    sys.stdout.flush()
forceMatchProcesses.postProcess(variables)

# terminate program clock
end = time.time()               # terminate program runtime timer
t = end - start                 # compute total running time
print "\nTotal running time: {0:.3f} sec".format(t)      # print total running time

# generate output files
outFile = open(variables.outFileName, 'w')
giveOutFiles(outFile)

# plot finalized data
forceMatchPlot.plotAllData(variables)