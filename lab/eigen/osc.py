#! /usr/bin/python

import sys
from math import sqrt,exp
from scipy import array,pi,factorial
from scipy.linalg import eig
from scipy.special import hermite

SIZE       = 100
ALPHA      = 0.1
HBAR_OMEGA = 1.0  # The same thing as E

def hamiltonian(size):
    a = array([[0. for i in xrange(size)] for i in xrange(size)])
    for i in xrange(size):
        for j in xrange(size):
            if i==j:
                a[i][j] += HBAR_OMEGA * (j+0.5) + \
                           -HBAR_OMEGA/4*(2*j+1) + \
                           ALPHA*HBAR_OMEGA/8*( 6*j**2 + 6*j + 3 )
            elif i==j+2:
                a[i][j] += -HBAR_OMEGA/4*sqrt(j+1)*sqrt(j+2) +\
                           ALPHA*HBAR_OMEGA/8*( sqrt(j+1)*sqrt(j+2)*(3+2*j) + \
                                     sqrt(j+1)*(j+2)**1.5+(j+1)**1.5*sqrt(j+2))
            elif i==j+4:
                a[i][j] += ALPHA*HBAR_OMEGA/8 * sqrt(j+1)*sqrt(j+2)*sqrt(j+3)*sqrt(j+4)
            elif i==j-2:
                a[i][j] += -HBAR_OMEGA/4*sqrt(j)*sqrt(j-1) +\
                           ALPHA*HBAR_OMEGA/8*( sqrt(j-1)*sqrt(j)*(2*j-1) + \
                                     sqrt(j-1)*j**1.5+(j-1)**1.5*sqrt(j) )
            elif i==j-4:
                a[i][j] += ALPHA*HBAR_OMEGA/8 * sqrt(j)*sqrt(j-1)*sqrt(j-2)*sqrt(j-3)
    return a

def sho_psi(n,xi):
    coeff = 1/sqrt( 2**n * factorial(n) * sqrt(pi) )
    expo  = exp(-1/2.* xi**2)
    herm  = hermite(float(n))(xi)
    return coeff*expo*herm

def new_psi(n,xi,E):
    ret = 0
    cutoff = 0
    for n,c in enumerate(E[n][1]): # the eigenvector
        cutoff += 1
        ret += c * sho_psi(n,xi)
        if cutoff==10:
            break
    return ret
        

def main():
    A = hamiltonian(SIZE)
    tmp = eig(A)

    def cmp( a, b ): #for eigenvalue/vector pairs
        return -1 if a[0]<b[0] else 0 if a[0]==b[0] else 1

    E = [0]*SIZE
    for i in xrange(SIZE):
        E[i] = ( tmp[0][i], [tmp[1][j][i] for j in xrange(SIZE)] )
    E.sort(cmp)

#     print "EIGENVALUES:"
#     for i in E:
#         print "%10.4f"% i[0]

    for i in range(40):
        sys.stderr.write("%i: " % i )
        print new_psi(0,-2+i/10.,E)
        #print sho_psi(3,-2+i/10.)
main()
