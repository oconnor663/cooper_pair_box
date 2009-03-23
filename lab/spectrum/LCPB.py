#! /usr/bin/python

import sys,os
from scipy import array, sqrt, exp, pi, factorial, cos, sin
from scipy.linalg import eig
from scipy.special import genlaguerre, poly1d

####
#### Create a function that gives an array of genlaguerres.
#### it will need EL and EC as parameters (for phi0), in
#### addition to a size
#### 

def ineff_test( phi0, size ):
    ret = [ [genlaguerre(i,j)(phi0**2/2) for j in range(size)]
            for i in range(size) ]
    return ret

def genlaguerre_array( phi0, size ):
    arg = phi0**2/2

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

def hamiltonian(size, EC, EJ, EL, flux):
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
                    * genlaguerre(n,2*m)(phi0**2/2)
            #the nonzero sine elements
            else:
                ### IS THIS PART RIGHT?
                n = min(row,col)
                m = (abs(col-row)-1)/2
                a[row][col] += \
                    -EJ * sin(2*pi*flux/flux0) * (-2)**(-m) * 2**-.5 \
                    * sqrt(factorial(n)/factorial(n+2*m+1)) \
                    * phi0**(2*m+1) * exp(phi0**2/-4) \
                    * genlaguerre(n,2*m+1)(phi0**2/2) ## Check overall signs
            
    return a

def sorted_eig( array ): ### real values only...
    vals, vecs = eig(array)
    vals = [ i.real for i in vals ]
    ret = zip( vals, vecs.T )
    def cmp( a, b ):
        return 0 if a[0]==b[0] else -1 if a[0]<b[0] else 1
    ret.sort(cmp=cmp)
    return ret

def main():

    EC = 1
    EJ = 1
    EL = 1
    FLUX = .25

    matrix_size = 50
    num_levels = 5
    
    print "Populating Matrix..."
    A = hamiltonian(matrix_size,EC,EJ,EL,FLUX)
    print "Diagonalizing..."
    print A
    e = sorted_eig(A)
    print "\nEnergies:"
    for i in e[:num_levels]: print i[0]
    

if __name__=='__main__':
    main()
