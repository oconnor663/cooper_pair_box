#! /usr/bin/python

import os,sys

def read_data( file ):
    fluxes = []
    data = []

    line = file.readline()
    while line:
        fluxes.append( float(line) )
        line = file.readline()
        data.append([])
        while line!='\n':
            data[-1].append(float(line))
            line = file.readline()
        line = file.readline()

    return (fluxes,data)

def main():
    (fluxes,data) = read_data( sys.stdin )

    out = os.popen("gnuplot","w")
    out.write( "set key off\n"
               "set terminal png\n"
               "set output 'out.png'\n"
               "plot '-' with points lt 3 pt 0\n" )
    for i in range(len(fluxes)):
        for j in data[i]:
            out.write( "%f %f\n" % (fluxes[i],j) )
    out.write( 'e\n' )
    out.close()

main()
