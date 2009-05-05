#! /usr/bin/python

import os,sys

def main():
    points = [ [float(i) for i in line.split()] for line in sys.stdin ]
    rows = len(points)
    cols = len(points[0])
    
    xmin, xmax, ymin, ymax = [float(i) for i in sys.argv[1:]]

    plotter = os.popen( "gnuplot", "w" )
    plotter.write( "set terminal png\n set output 'out.png'\n" )
    plotter.write( "set xrange [%f:%f]\n" % (xmin,xmax) )
    plotter.write( "set yrange [%f:%f]\n" % (ymin,ymax) )
    plotter.write( "set xlabel 'Magnetic flux'\n" )
    plotter.write( "set ylabel 'Energy (GHz)'\n" )
    plotter.write( "set key off\n set pm3d map\n" )
    plotter.write( "splot '-'\n" )

    point = 0

    for j in range(len(points[0])):
        y = ymin + (ymax-ymin)*j/cols
        for i in range(len(points)):
            point += 1
            x = xmin + (xmax-xmin)*i/rows
            if point>4250200:
                print ( "%i %f %f %f" % (point,y,x,points[i][j]) )
            plotter.write( "%f %f %f\n" % (y,x,points[i][j]) )
        plotter.write('\n')
    plotter.write('e\n')
    plotter.close()



if __name__=='__main__':
    main()

