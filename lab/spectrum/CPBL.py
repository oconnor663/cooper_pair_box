#! /usr/bin/python

import sys,os
from scipy import array, sqrt, exp, pi, factorial, cos, sin
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

def coshamiltonian( preham, flux ):
    flux0 = 1
    size = len(preham)
    ret = preham.copy()
    for row in range(size):
        for col in range(size):
            if (col-row)%2==0:
                ret[row][col] *= cos(2*pi*flux/flux0)
            else:
                ret[row][col] *= sin(2*pi*flux/flux0)
    return ret

def hamiltonian(cosham, EC, EJ, EL):
    hbar_w0 = sqrt( 8. * EL * EC )
    size = len(cosham)

    ret = cosham.copy()
    for row in range(size):
        for col in range(size):
            if (col-row)%2==0:
                ret[row][col] *= EJ
            else:
                ret[row][col] *= EJ
            if row==col:
                ret[row][col] += hbar_w0 * (row+0.5)            
    return ret

def sorted_eig( array ): ### real values only...
    vals, vecs = eig(array)
    vals = [ i.real for i in vals ]
    ret = zip( vals, vecs.T )
    def cmp( a, b ):
        return 0 if a[0]==b[0] else -1 if a[0]<b[0] else 1
    ret.sort(cmp=cmp)
    return ret

def make_curves( genlags, fluxes, EC, EJ, EL, num_curves=5 ):
    preham = prehamiltonian( genlags, EC, EJ, EL )
    ycurves = [ [ 0 for j in range(num_curves) ] for f in fluxes ]
    for i,flux in enumerate(fluxes):
        cosham = coshamiltonian( preham, flux )
        H = hamiltonian(cosham,EC,EJ,EL)
        e = sorted_eig(H)
        for j in range(num_curves):
            ycurves[i][j] = e[j+1][0]-e[0][0]

    return ycurves

def quad_diff( points, curves ):
    ### ASSUMES both generated over the same fluxes
    if len(points)!=len(curves):
        raise RuntimeError, "Something went wrong."
    
    sum = 0
    for i,yvals in enumerate(points):
        for data in yvals:
            sum += min( (data-theory)**2 for theory in curves[i] )
    return sum


calls = 0
def optimizer( EC_EJ_EL_tup, fluxes, points, genlags ):
    global calls
    EC, EJ, EL = EC_EJ_EL_tup 
    calls += 1
    curves = make_curves( genlags, fluxes, EC, EJ, EL )
    #plot_curves(fluxes,curves)
    ret = quad_diff( points, curves )
    print "Optimizer called #%i (val: %f)\n\t%.20f\n\t%.20f\n\t%.20f" \
         % (calls, ret, EC, EJ, EL)
    return ret

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

def guess_range(EC, EJ, EL):
    EC *= (.9+.2*random())
    EJ *= (.9+.2*random())
    EL *= (.9+.2*random())
    factor = 1.2  # This better be greater than SOMETHING b/c of preceding
    factor = float(factor)
    ret = ((EC/factor,factor*EC), 
           (EJ/factor, factor*EJ),
           (EL/factor, factor*EL))
    print "Guess ranges:"
    for i in ret: print i
    print "\n----------------------\n"
    return ret

def main():

    EC = 2.5
    EJ = 8.8
    EL = 0.5
    
    MATRIX_SIZE = 5
    NUM_POINTS = 10

    genlags = genlaguerre_array( MATRIX_SIZE )

    p = prehamiltonian(genlags,EC,EJ,EL)
    c = coshamiltonian(p,0)
    h = hamiltonian(c,EC,EJ,EL)
    print h

    fluxes = [ i*1./NUM_POINTS - .5 for i in range(NUM_POINTS) ]

    curves = make_curves( genlags, fluxes, EC, EJ, EL )
    
    #plot_curves( fluxes, curves )

    ranges = guess_range(EC,EJ,EL)

    print fmin_l_bfgs_b( optimizer, (0,0,0), None, (fluxes,curves,genlags),
                         True, ranges, factr=1e15)



if __name__=='__main__':
    main()
