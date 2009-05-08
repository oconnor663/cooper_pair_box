#! /usr/bin/python

import os,sys

def main():
    print "Reading data..."
    points = [ [float(i) for i in line.split()] for line in sys.stdin ]
    
    print "Paring data..."
    FAC = 1
    points = [ [points[i][j] for j in range(len(points[i])) if j%FAC==0]
               for i in range(len(points)) if i%FAC==0 ]
 
    rows = len(points)
    cols = len(points[0])

    
    ### Use this to input values for these parameters instead
    #xmin, xmax, ymin, ymax = [float(i) for i in sys.argv[1:]]
    xmin, xmax, ymin, ymax = -0.6414, 0.1151, 8.3, 11.799


    plotter = os.popen( "gnuplot", "w" )
    plotter.write( u"set xlabel 'Magnetic flux (\u03a6_0)'\n"
                   "set ylabel 'Energy (GHz)'\n".encode("utf-8"))
    plotter.write( "set terminal png enhanced font '/usr/share/fonts/truetype/freefont/FreeSerif.ttf' 14\n" )
    plotter.write( "set output 'out.png'\n" )
    plotter.write( "set xrange [%f:%f]\n" % (xmin,xmax) )
    plotter.write( "set yrange [%f:%f]\n" % (ymin,ymax) )
    plotter.write( "set key off\n set pm3d map\n" )
    plotter.write( "unset colorbox\n" )
    plotter.write( "set palette defined (-200 'green',-80 'purple', -40 'red',"
                   "-20 'orange', -10 'green', 0 'red', 10 'blue',"
                   "20 'purple', 40 'red', 70 'yellow',245 'green') \n" )
    plotter.write( "splot '-'\n" )

    print "Plotting..."

    for i in range(rows):
        # note that this is y, not x. To flip the image.
        y = ymin + (ymax-ymin)*i/rows
        for j in range(cols):
            x = xmin + (xmax-xmin)*j/cols
            plotter.write( "%f %f %f\n" % (x, y, points[i][j]))
        plotter.write('\n')
    plotter.write('e\n')
    plotter.close()



if __name__=='__main__':
    main()

