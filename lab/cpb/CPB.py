#! /usr/bin/python

import sys,os
from math import sqrt,exp
from scipy import array
from scipy.linalg import eig

MATRIX_RADIUS = 8

def hamiltonian(radius,ng,EJoverEC):
    size = 2*radius+1
    a = array([[0. for i in xrange(size)] for i in xrange(size)])
    for row in xrange(2*radius+1):
        for col in xrange(2*radius+1):
            if row==col:
                a[row][col] = 4.*(row-radius-ng)**2
            elif abs(row-col)==1:
                a[row][col] = -EJoverEC/2.
    return a

def make_plot( image_name, tmp_names, EJoverEC ):
    out = os.popen( "gnuplot", "w" )
    out.write( "set terminal png \n set output '%s'\n" % image_name )
    out.write( "set key off \n" )
    out.write( "set title 'Ej/Ec = %.2f' \n" % EJoverEC )
    out.write( "set yrange [-30:30] \n" )

    plot_command = "plot '%s' with lines" % tmp_names[0]
    for i in tmp_names[1:]:
        plot_command += ", '%s' with lines" % i
    out.write(plot_command)

    out.close()

def main():
    if len(sys.argv)<2:
        sys.stderr.write( "Supply a value for EJoverEC.\n" )
        sys.exit(1)
    EJoverEC = float(sys.argv[1])

    num_levels = 5
    tmp_names = [ os.popen("mktemp").read()[:-1] for i in range(num_levels) ]
    tmp_handles = [ open(i,"w") for i in tmp_names ]
    image_name = os.popen("mktemp").read()[:-1]
    os.system("mv %s %s.png" % (image_name,image_name))
    image_name += ".png"

    ng_rad = 7
    ng_pts = 1000
    for i in range(ng_pts):
        ng = -ng_rad + 2.*ng_rad/ng_pts*i
        A = hamiltonian(MATRIX_RADIUS,ng,EJoverEC)
        eigensystem = eig(A)
        #if i==500:
        #    print eigensystem
        energies = [ i.real for i in eigensystem[0] ]
        energies.sort()
        for j in range(num_levels):
            tmp_handles[j].write("%f %f\n" % (ng,energies[j]))
    for i in tmp_handles:
        i.close()

    make_plot(image_name,tmp_names,EJoverEC)
    os.system("kview %s" % image_name)


if __name__=='__main__':
    main()
