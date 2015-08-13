__author__ = 'Greg Poisson'

import numpy
import MDAnalysis

cF = None

def getFiles(vars):
    global cF
    cF = open(vars.configFile, 'r')
    print "\t'{}'".format(vars.configFile)

    getPsf(cF, vars)
    getForceDCDs(cF, vars)
    getCoordDCDs(cF, vars)
    getParam(cF, vars)
    getTemp(cF, vars)
    getGroups(cF, vars)

    initMDA(vars)

# Get name of PSF file from config file
def getPsf(cF, vars):
    while vars.psf is None:
        line = cF.next()
        if line == 'PSF FILE:\n':
            vars.psf = cF.next()[:-1]
            print('PSF File: {}'.format(vars.psf))
        elif line == 'END CONFIG FILE\n':
            print('No PSF file found in config.')
            break

# Get name of Force DCD files from config file
def getForceDCDs(cF, vars):
    vars.forceDcd = []
    txt = cF
    while len(vars.forceDcd) == 0:
        line = txt.next()
        if line == 'FORCE DCD FILES:\n':
            line = txt.next()
            while line != '\n':
                if line == 'END CONFIG\n':
                    print('NO FORCE DCD FILES FOUND')
                vars.forceDcd.append(line[:-1])
                line = txt.next()
            print('Force DCD files: {}'.format(vars.forceDcd))

# Get name of Coordinate DCD files from config file
def getCoordDCDs(cF, vars):
    vars.coordDcd = []
    txt = cF
    while len(vars.coordDcd) == 0:
        line = txt.next()
        if line == 'COORD DCD FILES:\n':
            line = txt.next()
            while line != '\n':
                if line == 'END CONFIG\n':
                    print('NO FORCE DCD FILES FOUND IN CONFIG')
                vars.coordDcd.append(line[:-1])
                line = txt.next()
            print('Coordinate DCD files: {}'.format(vars.coordDcd))

# Get name of the parameter file from config file
def getParam(cF, vars):
    txt = cF
    while vars.param is None:
        line = txt.next()
        if line == 'PARAMETER FILE:\n':
            line = txt.next()
            while line != '\n':
                if line == 'END CONFIG\n':
                    print('NO PARAMETER FILE FOUND IN CONFIG')
                vars.param = line[:-1]
                line = txt.next()
            print('Parameter file: {}'.format(vars.param))

# Set coordinate max/min and binsize
def getCoordBounds(cF, vars):
    txt = cF
    line = txt.next()
    length1 = len("MIN DISTANCE: ")
    length2 = len("MAX DISTANCE: ")
    length3 = len("BIN SIZE: ")

    # scan config file for coord and bin values
    while line != "END CONFIG\n":
        line = txt.next()
        if line[:length1] == "MIN DISTANCE: ":
            rem = -1 * (len(line) - length1)
            vars.rMin = int(line[rem:-1])
        elif line[:length2] == "MAX DISTANCE: ":
            rem = -1 * (len(line) - length2)
            vars.rMax = int(line[rem:-1])
        elif line[:length3] == "BIN SIZE: ":
            rem = -1 * (len(line) - length3)
            vars.binSize = float(line[rem:-1])

# Get temperature of system from config file
def getTemp(cF, vars):
    txt = cF
    while vars.temperature is None:
        line = txt.next().split()
        if len(line) > 1:
            if (line[0] == "SYSTEM") & (line[1] == "TEMPERATURE:"):
                vars.temperature = float(line[2])
                print "System temperature: {} K\n".format(vars.temperature)
            elif line[0] == "END":
                break

# Populate groups[] with atom groupings in config file
def getGroups(cF, vars):
    txt = cF
    while vars.groups is None:
        line = txt.next().split()
        if ("END" in line):
                break
        elif len(line) > 3:
            if (line[0] == "GIVE") & (line[1] == "ATOM") & (line[2] == "GROUPS:") & (line[3] == "ENABLED"):
                vars.grouping = True
                line = txt.next().split()
                vars.groups = []
                while len(line) > 0:
                    vars.group = []
                    for word in line:
                        # skip index number on each line
                        if ")" not in word:
                            # if values given are a range, add all included values
                            if "-" in word:
                                w = word.split("-")
                                w[1] = w[1].split(",")
                                begin = int(w[0])
                                end = int(w[1][0])
                                for i in range(begin, end + 1):
                                    vars.group.append(i)
                            # if value given has a comma on the end, remove it and
                            # append the value to the group
                            elif "," in word:
                                w = word.split(",")
                                atom = int(w[0])
                                vars.group.append(atom)
                            # append value to group
                            else:
                                vars.group.append(int(word))
                    vars.groups.append(vars.group)
                    line = txt.next().split()
            elif (line[0] == "GROUP") & (line[1] == "BY") & (line[2] == "RESIDUE:") & (line[3] == "ENABLED"):
                vars.groupingByResName = True
                line = txt.next().split()
                vars.groups = [line[3]]              # Set size of [groups] equal to config "Num of residues"

def initMDA(vars):
    vars.m = False       # Swtiched to True if user requests an rMax value greater than the system allows
    vars.force = MDAnalysis.Universe(vars.psf, vars.forceDcd)
    # coordinate universe
    vars.coord = MDAnalysis.Universe(vars.psf, vars.coordDcd)
    vars.dims = [vars.coord.dimensions[0], vars.coord.dimensions[1], vars.coord.dimensions[2]]
    vars.rMaxLimit = numpy.sqrt((vars.dims[0]**2) + (vars.dims[1]**2) + (vars.dims[2]**2))
    if vars.rMax > vars.rMaxLimit:
        vars.rMax = vars.rMaxLimit
        vars.m = True

    vars.binCount = int((vars.rMax - vars.rMin)/vars.binSize)