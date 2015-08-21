__author__ = 'Greg Poisson'

import numpy


# Define which plots to draw, and in which order
def getPlotOrder(cF, vars):
    txt = open(cF, 'r')
    line = txt.next().split()
    while len(line) == 0:
        line = txt.next().split()
    read = False
    while line[0] != "END":
        if (line[0] == "PLOTS") & (line[1] == "PLOT") & (line[2] == "ORDER"):   # Ignore config lines before plot order
            read = True
            line = txt.next().split()
        if read:                        # Begin reading plot order section of config
            count = int(line[-1:][0])   # Read user-defined number
            vars.plotOrder.append(count)
            line = txt.next().split()
            if len(line) == 0:
                break
        else:
            line = txt.next().split()
            while len(line) < 3:
                line = txt.next().split()

            if line[0] == "END":
                break

def executeProcesses(vars, firstFrame=False):
    # execute required processes
    for a in vars.atomGroups:
        for b in vars.atomGroups:
            if a.number < b.number:
                if firstFrame:
                    if pairIsUnique(a, b, vars):
                        vars.plots.append("{} {}".format(a.name, b.name))
                        vars.plots.append(buildBlankDataSet(a, b, vars))
                        computeCoordData(a, b, vars)
                        computeForceData(a, b, vars)
                    else:
                        computeCoordData(a, b, vars)
                        computeForceData(a, b, vars)
                else:
                    computeCoordData(a, b, vars)
                    computeForceData(a, b, vars)

# Identifies current particle pair, returns True if pair is unique (doesn't yet exist in plots array), False otherwise
def pairIsUnique(a, b, vars):
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in vars.plots:
        return False
    elif pairFlipped in vars.plots:
        return False
    else:
        return True

# Builds a blank data set in which to accumulate data
def buildBlankDataSet(a, b, vars):

    totForce = numpy.zeros(vars.binCount)
    integF = numpy.zeros(vars.binCount)

    freeEnergy = numpy.zeros(vars.binCount)
    rdf = numpy.zeros(vars.binCount)
    probDensity = numpy.zeros(vars.binCount) # Probability density of data counts
    dataCounts = numpy.zeros(vars.binCount)  # Keep track of the number of measurements made
                                             # in each bin, for averaging purposes

    dataSet = [totForce, integF, freeEnergy, rdf, probDensity, dataCounts]
    return dataSet

# Performs a series of coordinate-related computations on one pair of particles for one time step and returns the data
def computeCoordData(a, b, vars):
    magR = computeMagR(a, b, vars)[0]            # Determine distance between particles
    if vars.rMin <= magR < vars.rMax:             # If distance is within specified range
        binNo = int(magR / vars.binSize)     #   determine the appropriate bin number
        dataSet = findPair(a, b, vars)
        if binNo < vars.binCount:
            vars.plots[dataSet][len(vars.plots[dataSet]) - 1][binNo] += 1     # Increase measurement count by 1

# Performs a series of force computations on one pair of particles for one time step
def computeForceData(a, b, vars):
    magR = computeMagR(a, b, vars)[0]
    r = computeMagR(a, b, vars)[1]
    if vars.rMin < magR <= vars.rMax:
        binNo = int(magR / vars.binSize)
        if binNo < vars.binCount:
            dataSet = findPair(a, b, vars)
            vars.plots[dataSet][0][binNo] += computeMeanForce(magR, r, a, b)

# Determine magnitude of distance R between two particles
def computeMagR(a, b, vars):
    r = [0, 0, 0]
    for i in range(3):
        temp = a.position[i] - b.position[i]
        if temp < (-vars.dims[i])/2:
            temp += vars.dims[i]
        elif temp > (vars.dims[i]/2):
            temp -= vars.dims[i]
        r[i] = temp
    return numpy.linalg.norm(r), r

# Returns the index of the data set of a given pair of atoms, assuming it exists in the plots array
def findPair(a, b, vars):
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in vars.plots:
        return vars.plots.index(pair) + 1
    elif pairFlipped in vars.plots:
        return vars.plots.index(pairFlipped) + 1

# Perform post-data mining calculations
def postProcess(vars):
    print ("Done.\nConvert running sums to averages...")
    averageAll(vars)
    print ("Done.\nCompute RDF data...")
    rdf(vars)
    print ("Done.\nCompute Free Energy data...")
    freeEnergy(vars)
    print ("Done.\nCompute Potential of Mean Force data...")
    integrateForce(vars)
    print ("Done.\nCompute Distance Distribution data...")
    distanceDistribution(vars)
    print ("Done.\nConvert unmeasured values to NaN...")
    zeroToNan(vars)
    print ("Done.\n")

# Converts all running sums to averages
def averageAll(vars):
    lenPlots = len(vars.plots)
    for i in range(0, lenPlots):                  # Sort through all particle pairings
        if i % 2 == 1:
            for f in range(0, len(vars.plots[i])-1):
                for c in range(0, len(vars.plots[i][f])):   # Sort through all data sets
                    if vars.plots[i][-1:][0][c] != 0:
                        vars.plots[i][f][c] /= vars.plots[i][-1:][0][c]  # Divide running sum by measurement count

# Sets all uncomputed zeros to NaN
def zeroToNan(vars):
    plotsLength = len(vars.plots)
    setLength = len(vars.plots[1])
    subsetLength = len(vars.plots[1][0])
    for set in range(1, plotsLength):
        if set % 2 == 1:        # Skip over name elements, go into data elements
            for e in range(0, setLength - 1):       # Iterate through data sets
                for m in range(0, subsetLength - 1):       # Iterate through measurement counts
                    if vars.plots[set][setLength - 1][m] == 0:          # If there are zero measurements taken for a bin...
                        for q in range(0, setLength - 1):
                            vars.plots[set][q][m] = numpy.nan           # Set that bin = NaN for all subsets

# Integrates a set of force data for a given pair of atoms
def integrateForce(vars):
    for a in vars.atomGroups:
        for b in vars.atomGroups:
            if a.number < b.number:
                index = findPair(a, b, vars)
                sum = 0

                # Integrate the force data array, store in integrated force data array
                for tf in range(0, len(vars.plots[index][0]) - 1):
                    tf = len(vars.plots[index][0]) - 1 - tf
                    sum += vars.plots[index][0][tf] * vars.binSize
                    vars.plots[index][len(vars.plots[index])-5][tf] = sum

# Use force dcd file to examine total force interactions on particles
def computeMeanForce(magR, r, a, b):
    rHat = r / magR

    avgForce = (a.force - b.force) / 2
    MagAvgF = numpy.dot(avgForce, rHat)
    return MagAvgF

# Determine radial distribution frequency data
def rdf(vars):
    numSets = len(vars.plots)
    count = 0
    rho = getRho(vars)
    for set in range(0, numSets):
        if set % 2 == 1:
            for bin in range(1, vars.binCount):

                radius = bin * vars.binSize
                deltaR = vars.binSize
                volume = 4 * numpy.pi * radius**2 * deltaR

                hr = vars.plots[set][len(vars.plots[set])-1][bin]                     # Occurrences at distance r
                g = hr / (volume * rho * len(vars.coord.trajectory))             # g(r)
                vars.plots[set][len(vars.plots[set])-4][bin] = g                      # Send computed value to dataset
                count += 1


# Computes rho for RDF routine
def getRho(vars):
    ions = []
    for ion in vars.ionsCoord:
        if ion.name not in ions:
            count = 1
            ions.append(ion.name)
            ions.append(count)
        else:
            for i in range(len(ions)):
                if ion.name == ions[i]:
                    ions[i + 1] += 1
    N = 1
    for e in range(len(ions)):
        if e % 2 == 1:
            N *= ions[e]

    V = float(vars.dims[0]) * float(vars.dims[1]) * float(vars.dims[2])

    rho = N / V
    return rho

# Compute Free Energy data from Radial Distribution Data
def freeEnergy(vars):
    numSets = len(vars.plots)
    for set in range(0, numSets):
        if set % 2 == 1:
            for bin in range(0, vars.binCount):
                rdf = vars.plots[set][len(vars.plots[set])-4][bin]
                if rdf != 0.0:
                    logGr = numpy.log(rdf)         # Get probability and take log
                    fe = -vars.boltzmann * vars.temperature * logGr          # Compute free energy
                    vars.plots[set][len(vars.plots[set])-3][bin] = fe        # Send computed value to data set    `

# Convert distance frequency to distance probability distribution
def distanceDistribution(vars):
    lastSet = len(vars.plots)                             # index of last plots[] element (a set)
    countSubsetIndex = len(vars.plots[1])-1               # index of last subset in a set, containing
                                                          # observation counts at each bin length
    probDistIndex = len(vars.plots[1])-2                  # index of subset given for storing probability density data

    for set in range(0, lastSet):
        if set % 2 == 1:                             # set will reference all interaction pair datasets
            a = vars.plots[set][countSubsetIndex]
            b = numpy.sum(vars.plots[set][countSubsetIndex])
            n = a/b
            decimal = n / vars.binSize
            vars.plots[set][probDistIndex] = decimal