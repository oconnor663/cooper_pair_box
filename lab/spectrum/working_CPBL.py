#! /usr/bin/python

import sys,os
from scipy import array, sqrt, exp, pi, factorial, cos, sin
from scipy.linalg import eig
from scipy.special import genlaguerre, poly1d
from scipy.optimize import *
from random import random

def genlaguerre_array( arg, size ):
    ret = [ [ 0 for a in range(size) ] for n in range(size) ]
    for a in range(size):
        ret[0][a] = genlaguerre(0,a)
        ret[1][a] = genlaguerre(1,a)
        for n in range(2,size):
            ret[n][a] = ((2*(n-1)+a+1-poly1d([1,0]))*ret[n-1][a] - \
                             (n-1+a)*ret[n-2][a] ) /n
        for n in range(size):
            ret[n][a] = ret[n][a](arg)
    return ret

def hamiltonian(size, EC, EJ, EL, flux, genlags):
    hbar_w0 = sqrt( 8. * EL * EC )
    phi0 = ( 8. * EC / EL ) ** .25 
    flux0 = 1

    a = array([[0. for i in xrange(size)] for i in xrange(size)])
    for row in xrange(size):
        for col in xrange(size):
            # initialize diagonal elements
            if row==col:
                a[row][col] += hbar_w0 * (row+0.5)
            #the nonzero cosine elements
            if (col-row)%2==0:
                n = min(row,col)
                m = abs(col-row)/2 # because of Hermitianness
                a[row][col] += \
                    -EJ * cos(2*pi*flux/flux0) * (-2)**-m \
                    * sqrt(factorial(n)/factorial(n+2*m)) \
                    * phi0**(2*m) * exp(phi0**2/-4) \
                    * genlags[n][2*m]
            #the nonzero sine elements
            else:
                ### IS THIS PART RIGHT?
                n = min(row,col)
                m = (abs(col-row)-1)/2
                a[row][col] += \
                    -EJ * sin(2*pi*flux/flux0) * (-2)**(-m) * 2**-.5 \
                    * sqrt(factorial(n)/factorial(n+2*m+1)) \
                    * phi0**(2*m+1) * exp(phi0**2/-4) \
                    * genlags[n][2*m+1] ## Check overall signs
            
    return a

def sorted_eig( array ): ### real values only...
    vals, vecs = eig(array)
    vals = [ i.real for i in vals ]
    ret = zip( vals, vecs.T )
    def cmp( a, b ):
        return 0 if a[0]==b[0] else -1 if a[0]<b[0] else 1
    ret.sort(cmp=cmp)
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

def make_curves( fluxes, EC, EJ, EL, num_curves=5, matrix_size=20 ):
    phi0 = ( 8. * EC / EL ) ** .25 
    genlags = genlaguerre_array( phi0**2/2, matrix_size )
    
    ycurves = [ [ 0 for j in range(num_curves) ] for f in fluxes ]
    for i,flux in enumerate(fluxes):
        A = hamiltonian(matrix_size,EC,EJ,EL,flux,genlags)
        e = sorted_eig(A)
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
def optimizer( EC_EJ_EL_tup, fluxes, points ):
    global calls
    EC, EJ, EL = EC_EJ_EL_tup 
    calls += 1
    curves = make_curves( fluxes, EC, EJ, EL )
    #plot_curves(fluxes,curves)
    ret = quad_diff( points, curves )
    print "Optimizer called #%i (val: %f)\n\t%.20f\n\t%.20f\n\t%.20f" \
         % (calls, ret, EC, EJ, EL)
    return ret

def main():

    EC = 2.5
    EJ = 8.8
    EL = 0.5
    
    def guess_range(EC, EJ, EL):
        EC *= (.9+.2*random())
        EJ *= (.9+.2*random())
        EL *= (.9+.2*random())
        factor = 1.2  # This better be greater than SOMETHING b/c of preceding
        factor = float(factor)
        return ((EC/factor,factor*EC),
                (EJ/factor, factor*EJ),
                (EL/factor, factor*EL))

    num_points = 10
    fluxes = [ i*1./num_points - .5 for i in range(num_points) ]

    curves = make_curves( fluxes, EC, EJ, EL )
    
    #plot_curves( fluxes, curves )

    ranges = guess_range(EC,EJ,EL)
    print "Guess ranges:"
    for i in ranges: print i
    print "\n----------------------\n"

    print fmin_l_bfgs_b( optimizer, [0,0,0], None, (fluxes,curves),
                         True, ranges, factr=1e15)



#if __name__=='__main__':
#    main()

EC = 2.5
EJ = 8.8
EL = 0.5

phi0 = (8. * EC / EL)**.25
flux = .3345
H = hamiltonian( 20, EC, EJ, EL, flux, genlaguerre_array(phi0**2/2,20) )
e = sorted_eig(H)[:5]
print "\nEnergies:"
for i in [i[0] for i in e]:
    print i
    
