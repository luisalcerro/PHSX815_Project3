# Python code for 2D random walk. 
import sys
import numpy
import matplotlib.pyplot as plt
import pylab 
#from random import random
import math
#from scipy.special import factorial
sys.path.append(".")
from python.Random import Random




if __name__ == "__main__":
    if '-h' in sys.argv or '--help' in sys.argv:
        print ("Usage: %s [-seed number] [-Nwalk number] [-Nstep number] " % sys.argv[0])
        print
        sys.exit(1)

    # default seed
    seed = 786
    
    # default Nwalk
    Nwalk = 20000
    
    # default Nstep
    Nstep = 10

    # read the user-provided seed from the command line (if there)
    if '-seed' in sys.argv:
        p = sys.argv.index('-seed')
        seed = sys.argv[p+1]    
    if '-Nwalk' in sys.argv:
        p = sys.argv.index('-Nwalk')
        Nw = int(sys.argv[p+1])
        if Nw > 0:
            Nwalk = Nw
    if '-Nstep' in sys.argv:
        p = sys.argv.index('-Nstep')
        Ns = int(sys.argv[p+1])
        if Ns > 0:
            Nstep = Ns
    # class instance of our Random class using seed
    random = Random(seed)
    
    #############################################################
    ################ RANDOM WALK #######################
    #############################################################
    # Define positions and displacement in arrays
    x = numpy.zeros(Nstep) 
    y = numpy.zeros(Nstep)
    xfinal = []
    yfinal = []
    finald = []
    # filling the coordinates with Categorical numbers
    for j in range(1,Nwalk+1):
        for i in range(1, Nstep):
            z1, z2 = random.BoxMuller()
            x[i] = x[i - 1] + z1
            y[i] = y[i - 1] + z2
            finald.append(numpy.sqrt(x[Nstep-1]*x[Nstep-1]+y[Nstep-1]*y[Nstep-1]))
            xfinal.append(x[Nstep-1])
            yfinal.append(y[Nstep-1])
#########  PARAMETER ESTIMATION ############################
    addt = 0.
    for j in range(0,len(finald)):
        addt = addt + finald[j]*finald[j]
        par = numpy.sqrt(1./(2.*len(finald)) * addt)
################################################################
############### RANDOM WALK PLOTS#########################
################################################################
    plt.plot(xfinal, yfinal, 'p', color='steelblue')
    plt.grid(axis='x', alpha=0.75)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Final positions after ' + str(Nwalk) + ' walks and ' + str(Nstep) + ' steps')
    plt.show()
################################################################
    f1 = numpy.linspace(0.0,max(finald),1000)
    f2 = f1/(par**2) * numpy.exp(-f1*f1/(2.*par**2))
    plt.plot(f1,f2, color='black', linestyle='dashed', label = 'Rayleigh distribution, \u03C3= ' + str(round(par,4)))
    n, bins, patches = plt.hist(finald, 80, color ='steelblue', density = True, alpha=0.6, fill = True, hatch='/',histtype='step', linewidth=1, label = 'data')
    plt.legend(loc='upper right')
    plt.title('Distribution of the data compared to Rayleigh distribution ')
    plt.xlabel('r')
    plt.ylabel('Probability')
    plt.grid(axis='y', alpha=0.75)
    plt.show()
