## This is alex's test file

from optparse import OptionParser
import sys
from tools import *
from pysgpp import *
from painlesscg import cg,sd,cg_new
from math import sqrt
import random

from array import array

try:
    import psyco
    psyco.full()
    print "Using psyco"
except:
    pass
  
  
  

#-------------------------------------------------------------------------------
## Outputs a deprecated warning for an option
# @param option Parameter set by the OptionParser
# @param opt Parameter set by the OptionParser
# @param value Parameter set by the OptionParser
# @param parser Parameter set by the OptionParser
def callback_deprecated(option, opt, value, parser):
    print "Warning: Option %s is deprecated." % (option)
    

## Reads in an ARFF file:
# The data is stored in lists. There is a value list for every dimension of the data set. e.g. 
# [[2, 3],[1, 1]] are the data points P_1(2,1) and P_2(3,1)
#
# @param filename the file's filename that should be read
# @return returns a set of a array with the data (named data), a array with the classes (named classes) and the filename named as filename
def readDataARFF(filename):
    fin = open(filename, "r")
    data = []
    classes = []
    hasclass = False

    # get the different section of ARFF-File
    for line in fin:
        sline = line.strip().lower()
        if sline.startswith("%") or len(sline) == 0:
            continue

        if sline.startswith("@data"):
            break
        
        if sline.startswith("@attribute"):
            value = sline.split()
            if value[1].startswith("class"):
                hasclass = True
            else:
                data.append([])
    
    #read in the data stored in the ARFF file
    for line in fin:
        sline = line.strip()
        if sline.startswith("%") or len(sline) == 0:
            continue

        values = sline.split(",")
        if hasclass:
            classes.append(float(values[-1]))
            values = values[:-1]
        for i in xrange(len(values)):
            data[i].append(float(values[i]))
            
    # cleaning up and return
    fin.close()
    return {"data":data, "classes":classes, "filename":filename}

    
#-------------------------------------------------------------------------------
## Builds the training data vector
# 
# @param data a list of lists that contains the points a the training data set, coordinate-wise
# @return a instance of a DataVector that stores the training data
def buildTrainingVector(data):
    dim = len(data["data"])
    training = DataVector(len(data["data"][0]), dim)
    
    # i iterates over the data points, d over the dimension of one data point
    for i in xrange(len(data["data"][0])):
        for d in xrange(dim):
            training[i*dim + d] = data["data"][d][i]
    
    return training


def refinement2d():
    factory = Grid.createLinearBoundaryGrid(2)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(1)
    
    alpha = DataVector(9)
    alpha[0] = 0.0
    alpha[1] = 0.0
    alpha[2] = 0.0
    alpha[3] = 0.0
    alpha[4] = 0.0
    alpha[5] = 0.0
    alpha[6] = 0.0
    alpha[7] = 0.0
    alpha[8] = 1.0
    func = SurplusRefinementFunctor(alpha)
    
    #for n in xrange(factory.getStorage().size()):
    #    print factory.getStorage().get(n).getCoordinates()
    
    gen.refine(func)

    for n in xrange(factory.getStorage().size()):
        print factory.getStorage().get(n).getCoordinates()    


def refinement3d():
    factory = Grid.createLinearTrapezoidBoundaryGrid(3)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(1)
    
    alpha = DataVector(27)
    for i in xrange(len(alpha)):
         alpha[i] = 0.0

    alpha[26] = 1.0
    func = SurplusRefinementFunctor(alpha)
    
    gen.refine(func)
    
    for n in xrange(factory.getStorage().size()):
        print factory.getStorage().get(n).getCoordinates()  
   


def generateLaplaceMatrix(factory, level, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(level)
    
    laplace = factory.createOperationLaplace()
    
    # create vector
    alpha = DataVector(storage.size())
    erg = DataVector(storage.size())

    # create stiffness matrix
    m = DataVector(storage.size(), storage.size())
    m.setAll(0)
    for i in xrange(storage.size()):
        # apply unit vectors
        alpha.setAll(0)
        alpha[i] = 1
        laplace.mult(alpha, erg)
        if verbose:
            print erg, erg.sum()
        m.setColumn(i, erg)

    return m

        
def printDiagonal(m1):
    n = m1.getSize()

    # check diagonal
    values = []
    for i in range(n):
        values.append(m1[i*n + i])
    #values.sort()
    for i in range(n):
        print values[i]
        
        
def roundVector(v):
    for i in xrange(len(v)):
        v[i] = round(v[i], 10)
    
    return v

  
def printMatrix(m):
    n = m.getSize()

    for i in range(n):
        values = []
        for j in range(n):
            values.append(m[i*n + j])
        
        print roundVector(values)
        
    return
    
 
def checkSymmetry(m):
    error = 0.0
    tmp_error = 0.0
    n = m.getSize()

    for i in range(n):
        for j in range(n):
            tmp_error = 0.0
            tmp_error = m[i*n+j]-m[j*n+i]
            if abs(tmp_error) > abs(error):
                error = tmp_error
                
    return error
            
            
def printInnerPointsOnly(grid):
    pointstring = ""
    points = None
    printpoint = True
    
    for n in xrange(grid.getStorage().size()):
        pointstring = grid.getStorage().get(n).getCoordinates()
        points = pointstring.split()
        printpoint = True
        for i in xrange(len(points)):
            if points[i] == "0":
                printpoint = False
                
            if points[i] == "1":
                printpoint = False
                
        if (printpoint == True):
            print pointstring


def test_laplace():
    factory = Grid.createLinearBoundaryGrid(5)
    m = generateLaplaceMatrix(factory, 3)
    
    #print str(m) 
    #printDiagonal(m)
    
    print "C: Gitter mit Rand:"
    #printMatrix(m)
    print "\n"
    print "Symmetry-error:"
    print checkSymmetry(m)
    
    #factorytwo = Grid.createLinearGrid(2)
    #mtwo = generateLaplaceMatrix(factorytwo, 1)
    
    #print str(m) 
    #printDiagonal(m)
    
    #print "C down: Gitter ohne Rand:"
    #printMatrix(mtwo)
    
    return
                                  

def test_generation():
    factory = Grid.createLinearBoundaryGrid(3)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(2)
    
    #printInnerPointsOnly(factory)
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates()
        print points
    
    return

    
def test_refinement():
    factory = Grid.createLinearBoundaryGrid(2)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(0)
    
    alpha = DataVector(4)

    for i in xrange(len(alpha)):
        alpha[i] = 0.0

    alpha[0] = 1.0
    func = SurplusRefinementFunctor(alpha)
        
    gen.refine(func)
    
    alpha2 = DataVector(8)

    for i in xrange(len(alpha2)):
        alpha2[i] = 0.0

    alpha2[4] = 1.0
    func = SurplusRefinementFunctor(alpha2)
        
    gen.refine(func)
    
    alpha3 = DataVector(13)

    for i in xrange(len(alpha3)):
         alpha3[i] = 0.0

    alpha3[11] = 1.0
    func = SurplusRefinementFunctor(alpha3)
        
    gen.refine(func)
    
    alpha4 = DataVector(18)

    for i in xrange(len(alpha4)):
        alpha4[i] = 0.0

    alpha4[12] = 1.0
    func = SurplusRefinementFunctor(alpha4)
        
    gen.refine(func)
    
    alpha5 = DataVector(25)

    for i in xrange(len(alpha5)):
        alpha5[i] = 0.0

    alpha5[23] = 1.0
    func = SurplusRefinementFunctor(alpha5)    
    
    gen.refine(func)
        
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates()
        print points
    
    #print storage.size()
    
def test_new_refinement():
    factory = Grid.createLinearBoundaryGrid(3)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(0)
    
    alpha = DataVector(8)

    for i in xrange(len(alpha)):
        alpha[i] = 0.0

    alpha[0] = 1.0
    func = SurplusRefinementFunctor(alpha)
        
    gen.refine(func)
    
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates()
        print points    
    

def buildParableBoundary(dim):
    function = ""
    
    function = "1.0"
    
    for i in xrange(dim):
        function = function + "*(((-1.0)*(x" + str(i+1) + "-0.7)*(x" + str(i+1) + "-0.7))+2.0)"
        
    return function 


def buildLinearFunctionBoundary(dim):
    function = ""
    
    function = "1.0"
    
    for i in xrange(dim):
        function = function + "*(x" + str(i+1) + "+0.7)"
        
    return function 


def evalFunction(function, points):
    for i in xrange(len(points)):
        function = re.sub("x" + str(i+1), points[i], function)
            
    return eval(function)  


def doHierarchisation(node_values, grid):   
    tmp =  DataVector(grid.getStorage().size(), 1)
    
    for i in xrange(len(node_values)):
        tmp[i] = node_values[i]
    
    # create operation: hierarchisation
    hierarchisation = grid.createOperationHierarchisation()
    
    # execute hierarchisation
    hierarchisation.doHierarchisation(tmp)    

    return tmp

    
def test_getAffected():
    factory = Grid.createLinearBoundaryGrid(7)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(6)
 
    function = buildParableBoundary(7)
            
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates()
        print points 
        
    node_values = DataVector(storage.size(), 1)
    
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates().split()
        node_values[n] = evalFunction(function, points)
    #print node_values
    
    # do hierarchisation
    alpha = doHierarchisation(node_values, factory)
    
    p = DataVector(1,7)
    
    p[0] = 0.5
    p[1] = 0.5
    p[2] = 0.5
    p[3] = 0.5
    p[4] = 0.5
    p[5] = 0.5
    p[6] = 0.5
    pointstring = "0.5 0.5 0.5 0.5 0.5 0.5 0.5"
    evalpoints = pointstring.split()
    
    pc = factory.createOperationEval().eval(alpha, p)
    
    print "Value of Sparse Grid: " + str(pc)
    
    pc = evalFunction(function, evalpoints)
    
    print "Value of Function: " + str(pc)  
    
    
def test_grid_serialize():
    factory = Grid.createLinearBoundaryGrid(2)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(0)
    
    print factory.serialize()
    
    vector = DataVector(4)
    vector.setAll(0.0)
    vector[0] = 1.0
    func = SurplusRefinementFunctor(vector)
    gen.refine(func);
    
    print factory.serialize()


def openFile(filename):
    try:
        data = readDataARFF(filename)
    except:
        print ("An error occured while reading " + filename + "!")
        sys.exit(1)
        
    if data.has_key("classes") == False:
        print ("No classes found in the given File " + filename + "!")
        sys.exit(1)
        
    return data


def generateCMatrix(factory, level, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
    
    laplace = factory.createOperationLaplace()
    
    # create vector
    alpha = DataVector(storage.size())
    erg = DataVector(storage.size())

    # create stiffness matrix
    m = DataVector(storage.size(), storage.size())
    m.setAll(0)
    for i in xrange(storage.size()):
        # apply unit vectors
        alpha.setAll(0)
        alpha[i] = 1
        laplace.mult(alpha, erg)
        if verbose:
            print erg, erg.sum()
        m.setColumn(i, erg)

    return m


def generateBMatrix(factory, level, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    # create vector
    alpha = buildTrainingVector(openFile("function.out"))
    point = DataVector(storage.size())
    erg = DataVector(storage.size())

    # create stiffness matrix
    m = DataVector(storage.size(), storage.size())
    m.setAll(0)
    for i in xrange(storage.size()):
        # apply unit vectors
        point.setAll(0)
        point[i] = 1
        b.mult(alpha, point, erg)
        if verbose:
            print erg, erg.sum()
        m.setColumn(i, erg)

    return m


def generateBtransMatrix(factory, level, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    # create vector
    alpha = buildTrainingVector(openFile("function.out"))
    point = DataVector(storage.size())
    erg = DataVector(storage.size())

    # create stiffness matrix
    m = DataVector(storage.size(), storage.size())
    m.setAll(0)
    for i in xrange(storage.size()):
        # apply unit vectors
        point.setAll(0)
        point[i] = 1
        b.multTranspose(alpha, point, erg)
        if verbose:
            print erg, erg.sum()
        m.setColumn(i, erg)

    return m


def build_DM_Matrices():
    factory = Grid.createLinearBoundaryGrid(3)
    level = 3
    gen = factory.createGridGenerator()
    gen.regular(level)
    
    laplace_m = generateCMatrix(factory, level)
    B_m = generateBMatrix(factory, level)
    B_trans_m = generateBtransMatrix(factory, level)    
    

def calcNumberGridPoints():
    factory = Grid.createLinearBoundaryGrid(2)
    level = 6
    gen = factory.createGridGenerator()
    gen.regular(level)

    print factory.getStorage().size()    
    

#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    test_grid_serialize()
    #test_generation()