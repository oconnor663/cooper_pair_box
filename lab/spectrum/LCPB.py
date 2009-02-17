#! /usr/bin/python

import sys,os
from scipy import array, sqrt, exp, pi, factorial
from scipy.linalg import eig
from scipy.special import genlaguerre

def hamiltonian(size, EC, EJ, EL):
    hbar_w0 = sqrt( 8. * EL * EC )
    phi0 = ( 8. * EC / EL ) ** .25 

    a = array([[0. for i in xrange(size)] for i in xrange(size)])
    for row in xrange(size):
        for col in xrange(size):
            # initialize diagonal elements
            if row==col:
                a[row][col] += hbar_w0 * (row+0.5)
            #the nonzero cosine elements
            if (col-row)%2==0:
                n = row if row<col else col
                m = abs(col-row)/2 # because of Hermitianness
                a[row][col] += -EJ * (-2)**-m * sqrt(factorial(n)/factorial(n+2*m)) \
                               * phi0**(2*m) * exp(phi0**2/-4) \
                               * genlaguerre(n,2*m)(phi0**2/2)
            #the nonzero sine elements
            else:
                continue  # flux term needs to be fixed
                n = row # THIS IS WRONG
                m = (abs(col-row)-1)/2
                a[row][col] += -EJ * (-2)**-m*sqrt(factorial(n)/factorial(n+2*m+1)) \
                               * phi0**(2*m+1) * exp(phi0**2/-4) \
                               * genlaguerre(n,2*m+1)(phi0**2/2)
            
    print a
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
#     if len(sys.argv)<2:
#         sys.stderr.write( "Supply a value for EJoverEC.\n" )
#         sys.exit(1)
#     EJoverEC = float(sys.argv[1])

    EC = 1
    EJ = 1
    EL = 1

    num_levels = 5

#    for i in range(num_levels,20):
    for i in range(5,6):
        matrix_size = i
        A = hamiltonian(matrix_size,EC,EJ,EL)
        e = sorted_eig(A)
        print
        print i
        print
        for i in e[:num_levels]: print i[0]
        print
    

if __name__=='__main__':
    main()
