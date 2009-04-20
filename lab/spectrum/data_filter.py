#! /usr/bin/python

import sys,os
from math import sqrt

def outliers( col, tolerance ):
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
        if abs(col[i]-avg)>tolerance*dev:
            ret.append(i)
    return ret

def prune( points, min_neighbors ):

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

    sys.stderr.write( "Parsing...\n" )
    
    raw_data = [ [float(i) for i in line.split()]
                 for line in sys.stdin ]
    
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

    MIN_FLUX = -1.282
    MAX_FLUX = 0.224
    MIN_E = 8.5818
    MAX_E = 9.0548

    DIFF_FLUX = MAX_FLUX - MIN_FLUX
    DIFF_E = MAX_E - MIN_E

    ### Print euclidean indexing
    sys.stderr.write( "Printing...\n" )
    width = len(points)
    height = len(points[0])
    
    for i in range(width):
        #print flux
        print MIN_FLUX + DIFF_FLUX * i / width
        for j in range(height):
            #print energy values
            if points[i][j]==1 or points_copy[i][j]==1:
                print MIN_E + DIFF_E * j / height
        print

if __name__=="__main__":
    main()
