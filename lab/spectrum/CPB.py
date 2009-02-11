#! /usr/bin/python

import sys,os
from scipy import array, sqrt, exp, pi
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

def make_plot( image_name, tmp_names, overlay_names, EJoverEC ):
    out = os.popen( "gnuplot", "w" )
    out.write( "set terminal png \n set output '%s'\n" % image_name )
    out.write( "set key off \n" )
    out.write( "set title 'Ej/Ec = %.2f' \n" % EJoverEC )
    out.write( "set yrange [-30:30] \n" )
    
    FIRST_STYLE = "with lines lt 1 lw 2"
    SECOND_STYLE = "with lines lt 2 lw 1"
    plot_command = "plot '%s' %s" % (tmp_names[0],FIRST_STYLE)
    for i in tmp_names[1:]:
        plot_command += ", '%s' %s" % (i,FIRST_STYLE)
    for i in overlay_names:
        plot_command += ", '%s' %s" % (i,SECOND_STYLE)
    out.write(plot_command)

    out.close()

def sorted_eig( array ): ### real values only...
    vals, vecs = eig(array)
    vals = [ i.real for i in vals ]
    ret = zip( vals, vecs.T )
    def cmp( a, b ):
        return 0 if a[0]==b[0] else -1 if a[0]<b[0] else 1
    ret.sort(cmp=cmp)
    return ret

def phi_waveform( n_array, phi ):
    n_indices = range(-(len(n_array)/2), len(n_array)/2+1)
    if len(n_array)!=len(n_indices):
        raise RuntimeError, "Uh oh"
    ret = 0
    for n,coeff in zip(n_indices, n_array):
        ret += coeff * 1/sqrt(2*pi) * exp(1j*phi*n)
    return abs(ret)**2

def main():
    if len(sys.argv)<2:
        sys.stderr.write( "Supply a value for EJoverEC.\n" )
        sys.exit(1)
    EJoverEC = float(sys.argv[1])

    num_levels = 5

    tmp_names_en = [ os.popen("mktemp").read()[:-1] for i in range(num_levels) ]
    tmp_handles_en = [ open(i,"w") for i in tmp_names_en ]

    tmp_names_wave = [ os.popen("mktemp").read()[:-1] for i in range(num_levels) ]
    tmp_handles_wave = [ open(i,"w") for i in tmp_names_wave ]

    image_name = os.popen("mktemp").read()[:-1]
    os.system("mv %s %s.png" % (image_name,image_name))
    image_name += ".png"

    ng_rad = 2
    num_pts = 100
    for i in range(num_pts):
        ng = -ng_rad + 2.*ng_rad/num_pts*i
        A = hamiltonian(MATRIX_RADIUS,ng,EJoverEC)
        e = sorted_eig(A)
        for j in range(num_levels):
            tmp_handles_en[j].write("%f %f\n" % (ng,e[j][0]))
    for i in tmp_handles_en:
        i.close()

    phi_rad_factor = pi
    waveform_scale_factor = 5
    phi_ng_value = 0
    A = hamiltonian(MATRIX_RADIUS,phi_ng_value,EJoverEC)
    e = sorted_eig(A)

    for i in range(num_pts):
        phi = -ng_rad + 2.*ng_rad/num_pts*i # ng_rad used for axis consistency
        for j in range(num_levels):
            tmp_handles_wave[j].write("%f %f\n" % (phi,
                               waveform_scale_factor * \
                               phi_waveform(e[j][1],phi*phi_rad_factor) + \
                               e[j][0]) )
    for i in tmp_handles_wave:
        i.close()

    make_plot(image_name,tmp_names_en,tmp_names_wave,EJoverEC)
    os.system("display %s" % image_name)


if __name__=='__main__':
    main()
