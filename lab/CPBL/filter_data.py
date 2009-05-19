#! /usr/bin/python

import sys,os
from math import sqrt

def outliers( col, tolerance ):
    """Returns the points that are outside the average phase shift (z
 axis) range for a given column"""
    length = len(col)
    avg = 0
    avg2 = 0
    for i in col:
        avg+=i
        avg2+=i**2
    avg /= length
    avg2 /= length
    dev = sqrt(avg2 - avg**2)
    
    ret = []
    for i in range(length):
        #
        #  This is the formula for inclusion.
        #        
        if abs(col[i]-avg)>tolerance*dev:
            ret.append(i)
    return ret

def prune( points, min_neighbors ):
    """Returns points that have at least a minimum number of adjacent neighbors
    whoe also passed the previous test"""

    def good_point( point, height, width ):
        return point[0]>=0 and point[0]<height \
            and point[1]>=0 and point[1]<width
    
    HEIGHT = len(points)
    WIDTH = len(points[0])
    toprune = []
    for i in range(HEIGHT):
        for j in range(WIDTH):
            if points[i][j]==1:
                neighbors = 0
                for xoff in range(-1,2):
                    for yoff in range(-1,2):
                        if good_point( (i+xoff,j+yoff), HEIGHT, WIDTH ):
                            neighbors += points[i+xoff][j+yoff]
                if neighbors-1 < min_neighbors:
                    toprune.append((i,j))
    for i,j in toprune:
        points[i][j] = 0

def streak( points, min_streak ):
    """Remove points from any column unless they lie in a row at least
    min_streak points long in that column"""

    for col in points:
        L = len(col)
        a = 0
        b = 0
        while True:
            if b==L or col[b]==0:
                if b-a < min_streak:
                    for i in range(a,b):
                        col[i] = 0
            if b==L:
                break
            if col[b]==0:
                b += 1
                a = b
            else:
                b += 1

def main():

    TOLERANCE = 1.5
    MIN_NEIGHBORS = 6
    MIN_STREAK = 6

    if len(sys.argv)<4:
        sys.stderr.write( "Provide a data file, an energies file (for the y axis), and an outfile name as arguments.\n" )
        return 1

    datafile = open(sys.argv[1])
    yrangefile = open(sys.argv[2])
    outfile = open(sys.argv[3],'w')

    sys.stderr.write( "Parsing...\n" )
    
    raw_data = [ [float(i) for i in line.split()]
                 for line in datafile ]
    freqencies = [ float(i)/1e9 for i in yrangefile ]
    
    # rotate plot from matrix indexing to euclidean indexing
    sys.stderr.write( "Rotating...\n" )
    data = [ [ j[i] for j in raw_data ]
             for i in range(len(raw_data[0])) ]
    del raw_data

    ### Selecting only those points with value outside TOLERANCE

    sys.stderr.write( "Scanning...\n" )
    points = [ [0 for i in col] for col in data ]
    points_copy = [ [0 for i in col] for col in data ]
    for i,col in enumerate(data):
        for j in outliers(col,TOLERANCE):
            points[i][j] = 1
            points_copy[i][j] = 1

    ### Pruning isolated points
    sys.stderr.write( "Pruning...\n" )
    prune(points,MIN_NEIGHBORS)
    
    ### Streaking the copy
    sys.stderr.write( "Streaking...\n" )
    streak(points_copy,MIN_STREAK)

    #
    # NB: THESE NEED TO BE CHANGED FOR NEW DATA SETS
    #
    MIN_FLUX = -0.6414
    MAX_FLUX = 0.1151
    DIFF_FLUX = MAX_FLUX - MIN_FLUX

    ### Print euclidean indexing
    sys.stderr.write( "Printing...\n" )
    width = len(points)
    height = len(points[0])
    
    for i in range(width):
        #print flux
        outfile.write( "%f\n" % (MIN_FLUX + DIFF_FLUX * i / width) )
        for j in range(height):
            #print energy (frequency) values
            if points[i][j]==1 or points_copy[i][j]==1:
                outfile.write( "%f\n" % freqencies[j] )
        outfile.write("\n")

if __name__=="__main__":
    main()
