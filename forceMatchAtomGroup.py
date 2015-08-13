__author__ = 'Greg Poisson'

class fmAtomGroup():
    name = None
    number = 0
    position = [0, 0, 0]
    mass = 0.0
    charge = 0.0
    force = [0, 0, 0]

    # Call toString on an atom group for identifying information
    def toString(self):
        print "\tfmAtomGroup -- Name: {}\tNumber: {}".format(self.name, self.number)

# Select and perform a grouping operation during first timestep
def groupAtoms(vars):
    parseSolvent(vars)

    if vars.grouping:
        groupByCustom(vars)
        print "\tAtoms successfully grouped.\n\t\t{} atoms -> {} groups".format(len(vars.ionsCoord), len(vars.atomGroups))
    elif vars.groupingByResName:
        groupByResName(vars)
        print "\tAtoms successfully grouped.\n\t\t{} atoms -> {} groups".format(len(vars.ionsCoord), len(vars.atomGroups))
    else:
        noGrouping(vars)
        print "\tNo grouping was applied.\n\t\t{} atoms.".format(len(vars.atomGroups))


# Group atoms according to user-defined groups
def groupByCustom(vars):
    g = []
    print "\tCustom atom grouping requested..."
    count = 1
    for group in vars.groups:
        particle = fmAtomGroup()
        particle.name = str(count)
        particle.number = count
        particle.position = getGroupCOM(vars, group)
        particle.mass = getGroupMass(vars, group)
        particle.charge = getGroupCharge(vars, group)
        particle.force = getGroupForce(vars, group)
        count += 1

        g.append(particle)

    vars.atomGroups = g

# Group atoms according to their MDAnalysis resname
def groupByResName(vars):
    g = []
    print "\tAtom grouping by residue name requested..."
    rParticle = []
    currResName = vars.ionsCoord[0].resname
    for ion in vars.ionsCoord:
        if ion.resname == currResName:
            rParticle.append(ion.number)
        else:
            fParticle = fmAtomGroup
            fParticle.name = str(ion.resname)
            fParticle.number = len(g)
            fParticle.position = getGroupCOM(vars, rParticle)
            fParticle.mass = getGroupMass(vars, rParticle)
            fParticle.charge = getGroupCharge(vars, rParticle)
            fParticle.force = getGroupForce(vars, rParticle)

            g.append(fParticle)
            rParticle = []
            rParticle.append(ion.number)
            currResName = ion.resname

    vars.atomGroups = g

# Direct conversion MDAnalysis AtomGroup to fmAtomGroup
def noGrouping(vars):
    g =[]
    print "\tNo atom grouping requested."
    ionsCount = len(vars.ionsCoord)
    for ion in range(0, ionsCount):
        particle = fmAtomGroup
        particle.name = vars.ionsCoord[ion].name
        particle.number = ion
        particle.position = vars.ionsCoord[ion].position
        particle.mass = vars.ionsCoord[ion].mass
        particle.charge = vars.ionsCoord[ion].charge
        particle.force = vars.ionsForce[ion].position

        g.append(particle)

    vars.atomGroups = g

# Update each atomGroup's data
def updateGroups(vars):
    count = 0
    for g in vars.groups:
        vars.atomGroups[count].position = getGroupCOM(vars, g)
        vars.atomGroups[count].mass = getGroupMass(vars, g)
        vars.atomGroups[count].charge = getGroupCharge(vars, g)
        vars.atomGroups[count].force = getGroupForce(vars, g)
        count += 1

# Define subset of data without solvent
def parseSolvent(vars):
    vars.ionsCoord = vars.coord.selectAtoms("not name H1 and not name H2 and not name OH2")
    vars.ionsForce = vars.force.selectAtoms("not name H1 and not name H2 and not name OH2")

# Return center of mass of a group of atoms
def getGroupCOM(vars, group):
    sumOfMasses = 0
    sumProdcuts = 0
    count = 0
    for number in group:
        mi = vars.ionsCoord[count].mass
        ri = vars.ionsCoord[count].position
        sumOfMasses += mi
        sumProdcuts += (mi * ri)
        count += 1
    return ((1/sumOfMasses) * sumProdcuts)

# Return total mass of a group of atoms
def getGroupMass(vars, group):
    mass = 0
    count = 0
    for number in group:
        mass += vars.ionsCoord[count].mass
        count += 1
    return mass

def getGroupCharge(vars, group):
    charge = 0
    count = 0
    for number in group:
        charge += vars.ionsCoord[count].charge
        count += 1
    return float(charge / count)

def getGroupForce(vars, group):
    force = [0, 0, 0]
    count = 0
    for number in group:
        force += vars.ionsForce[count].position
        count += 1
    return force

# Checks to see if a string represents a float
def isFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False