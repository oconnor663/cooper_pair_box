#! /usr/bin/python

import sys,os
from math import sqrt
from scipy import array, exp, pi, factorial, cos, sin, dot
from scipy.linalg import eig
from scipy.special import genlaguerre, poly1d
from scipy.optimize import *
from random import random

def genlaguerre_array( size ):
    ret = [ [ 0 for a in range(size) ] for n in range(size) ]
    for a in range(size):
        ret[0][a] = genlaguerre(0,a)
        ret[1][a] = genlaguerre(1,a)
        for n in range(2,size):
            ret[n][a] = ((2*(n-1)+a+1-poly1d([1,0]))*ret[n-1][a] - \
                             (n-1+a)*ret[n-2][a] ) /n
    return ret

def prehamiltonian( genlags, EC, EJ, EL ):
    ### NB: omits a factor of EJ (on top of the flux dependency and the 
    ### diagonal terms) for use as a derivative later
    size = len(genlags)
    hbar_w0 = sqrt( 8. * EL * EC )
    phi0 = ( 8. * EC / EL ) ** .25 
    arg = phi0**2/2
    
    genlags = [[ f(arg) for f in row ] for row in genlags]
    ret = [ range(size) for i in range(size) ] #values set below

    for row in range(size):
        for col in range(size):
            #the nonzero cosine elements
            if (col-row)%2==0:
                n = min(row,col)
                m = abs(col-row)/2 # because of Hermitianness
                ret[row][col] = -(-2)**-m \
                    * sqrt(factorial(n)/factorial(n+2*m)) \
                    * phi0**(2*m) * exp(phi0**2/-4) \
                    * genlags[n][2*m]
            #the nonzero sine elements
            else:
                ### IS THIS PART RIGHT?
                n = min(row,col)
                m = (abs(col-row)-1)/2
                ret[row][col] = -(-2)**(-m) * 2**-.5 \
                    * sqrt(factorial(n)/factorial(n+2*m+1)) \
                    * phi0**(2*m+1) * exp(phi0**2/-4) \
                    * genlags[n][2*m+1] ## Check overall signs
    return array(ret)

def sorted_eig( array ): ### real values only...
    vals, vecs = eig(array)
    vals = [ i.real for i in vals ]
    ret = zip( vals, vecs.T )
    def cmp( a, b ):
        return 0 if a[0]==b[0] else -1 if a[0]<b[0] else 1
    ret.sort(cmp=cmp)
    return ret

SQRTS = [sqrt(i) for i in range(20)]

def opC( vec, negsqrtELo2EC ):
    size = len(vec)
    ret = array([0. for i in vec])
    for i,a in enumerate(vec):
        if i>1:
            ret[i-2] += SQRTS[i]*SQRTS[i-1]*a
        if i<size-2:
            ret[i+2] += SQRTS[i+1]*SQRTS[i+2]*a
        ret[i] += -(2*i+1)*a
    ret *= negsqrtELo2EC
    return ret

def opL( vec, sqrtECo2EL ):
    size = len(vec)
    ret = array([0. for i in vec])
    for i,a in enumerate(vec):
        if i>1:
            ret[i-2] += SQRTS[i]*SQRTS[i-1]*a
        if i<size-2:
            ret[i+2] += SQRTS[i+1]*SQRTS[i+2]*a
        ret[i] += (2*i+1)*a #here
    ret *= sqrtECo2EL #here
    return ret

def solve_energies( preham, EC, EJ, EL, flux, num ):
    '''Returns the energies and their derivatives for a single flux
    point, as a nested tuple. CALCULATES DIFFERENCE FROM GROUND.'''

    flux0 = 1
    hbar_w0 = sqrt( 8. * EL * EC )
    size = len(preham)
    if num+1>len(preham): raise RuntimeError, "mistake"
    cosham = preham.copy()
    for row in range(size):
        for col in range(size):
            if (col-row)%2==0:
                cosham[row][col] *= cos(2*pi*flux/flux0)
            else:
                cosham[row][col] *= sin(2*pi*flux/flux0)
    H = cosham.copy()
    for row in range(size):
        for col in range(size):
            if (col-row)%2==0:
                H[row][col] *= EJ
            else:
                H[row][col] *= EJ
            if row==col:
                H[row][col] += hbar_w0 * (row+0.5)
    system = sorted_eig(H)[:num+1]

    results = []
    negsqrtELo2EC = -sqrt(EL/(2.*EC))
    sqrtECo2EL = sqrt(EC/(2.*EL))
    for level in system:
        dEC = dot( level[1], opC(level[1],negsqrtELo2EC) )
        dEL = dot( level[1], opL(level[1],sqrtECo2EL) )
        dEJ = dot( level[1], dot(cosham,level[1]) )
        results.append( array([level[0], dEC, dEJ, dEL]) )
    return [ i-results[0] for i in results[1:] ]

queries = 0
def niceness( (EC,EJ,EL), genlags, fluxes, data, num_curves ):
    '''Returns the value of the difference of squares and the
    vector of partials in a format compatible with the optimizer.
    ( f, array(f_EC,f_EJ,f_EL) )'''

    global queries
    queries += 1
    print "Call #%i" %queries
    
    f = 0        # the sum of squares
    f_EC = 0     # and its partial derivatives
    f_EJ = 0
    f_EL = 0

    num_points = 0

    P = prehamiltonian(genlags,EC,EJ,EL)

    for i,flux in enumerate(fluxes):
        E = solve_energies(P,EC,EJ,EL,flux, num_curves )
        
        for d in data[i]:
            num_points += 1
            index = 0
            diffsq = (E[0][0] - d)**2
            # now find the energy it's closest to
            for j in range(1,num_curves):
                newdiff = (E[j][0] - d)**2
                if newdiff < diffsq:
                    diffsq = newdiff
                    index = j
            # now add the results to the running totals
            diff = E[index][0] - d
            f += diffsq
            f_EC += 2 * diff * E[index][1]
            f_EJ += 2 * diff * E[index][2]
            f_EL += 2 * diff * E[index][3]
    
    print "Niceness queried (# %i): %f" \
        "\n\t%.20f\n\t%.20f\n\t%.20f\n" % \
        (queries,f,EC,EJ,EL)

    f /= num_points
    f_EC /= num_points
    f_EJ /= num_points
    f_EL /= num_points

    return (f, array((f_EC,f_EJ,f_EL)))

def make_data( genlags, fluxes, EC, EJ, EL, num_curves ):
    P = prehamiltonian( genlags, EC, EJ, EL )
    data = [ [ 0 for j in range(num_curves) ] for f in fluxes ]
    for i,flux in enumerate(fluxes):
        E = solve_energies( P, EC, EJ, EL, flux, num_curves )
        for j in range(num_curves):
            data[i][j] = E[j][0]
    return data

def plot_curves( xpoints, ycurves ): # ycurves is a nested list
    # at the moment, ycurves contains sorted lists of each yval
    # at a given flux (i.e. not lists representing continuous curves)
    picname = os.popen( "mktemp", "r" ).read()[:-1]
    grapher = os.popen( "gnuplot", "w" )
    grapher.write( "set key off\n" )
    grapher.write( "set terminal png\nset output '%s'\n" % picname )
    grapher.write( "plot '-' w l" )
    for i in range(len(ycurves[0])-1):
        grapher.write( ", '-' w l" )
    grapher.write('\n' )
    for curve in range(len(ycurves[0])):
        for i,x in enumerate(xpoints):
            grapher.write( "%f %f\n" % (x,ycurves[i][curve]) )
        grapher.write('e\n')
    grapher.close()

    os.system( "display %s" % picname )

def guess_bounds(EC, EJ, EL):
    ret = ((EC*random(), (1+random())*EC), 
           (EJ*random(), (1+random())*EJ),
           (EL*random(), (1+random())*EL))
    return ret

def main():

    EC = 10*random() #2.5
    EJ = 10*random() #8.8
    EL = 10*random() #0.5

    print "EC = %f\nEJ = %f\nEL = %f" % (EC,EJ,EL)
    
    MATRIX_SIZE = 20
    NUM_ENERGIES = 5
    NUM_POINTS = 30

    genlags = genlaguerre_array( MATRIX_SIZE )
    fluxes = [ i*1./NUM_POINTS - .5 for i in range(NUM_POINTS) ]

    data = make_data( genlags,fluxes,EC,EJ,EL,NUM_ENERGIES )

    #plot_curves( fluxes, data )

    bounds = guess_bounds(EC,EJ,EL)
    guess = array( ( (bounds[0][0]+bounds[0][1])/2,
                     (bounds[1][0]+bounds[1][1])/2,
                     (bounds[2][0]+bounds[2][1])/2 ) )

    print "Guess ranges:"
    for i in bounds: print " ",i
    print "Initial guess:",guess
    
    ret = fmin_l_bfgs_b( niceness, guess,
                         args = (genlags,fluxes,data,NUM_ENERGIES),
                         bounds = bounds )

    print "\nNumber of calls:", queries # a global var
    print "EC = %f\nEJ = %f\nEL = %f" % (ret[0][0],ret[0][1],ret[0][2])

    if ret[1] > .01:
        raise RuntimeError, "Failure to converge."


if __name__=='__main__':
    main()


#g = genlaguerre_array(20)
#p = prehamiltonian( g, 2.5,8.8,0.5)
#solve_energies(p,2.5,8.8,0.5,.3345,5)
