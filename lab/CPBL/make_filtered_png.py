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

    ### Use this to input values for these parameters instead
    #xmin, xmax, ymin, ymax = [float(i) for i in sys.argv[1:]]
    xmin, xmax, ymin, ymax = -0.6414, 0.1151, 8.3, 11.799

    out = os.popen("gnuplot","w")
    out.write( u"set xlabel 'Magnetic flux (\u03a6_0)'\n"
               "set ylabel 'Energy (GHz)'\n".encode("utf-8"))
    out.write( "set key off\n"
               "set terminal png enhanced font '/usr/share/fonts/truetype/freefont/FreeSeriv.ttf' 14\n"
               "set output 'out.png'\n"
               "set xrange [%f:%f]\n"
               "set yrange [%f:%f]\n"
               "plot '-' with points lt 3 pt 0\n" %
               (xmin,xmax,ymin,ymax))
    for i in range(len(fluxes)):
        for j in data[i]:
            out.write( "%f %f\n" % (fluxes[i],j) )
    out.write( 'e\n' )
    out.close()

main()
