#! /usr/bin/python

import sys,copy
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

    if len(sys.argv)>1:
        TOLERANCE = float(sys.argv[1])
    else:
        TOLERANCE = 1.0

    if len(sys.argv)>2:
        MIN_NEIGHBORS = int(sys.argv[2])
    else:
        MIN_NEIGHBORS = 8

    if len(sys.argv)>3:
        MIN_STREAK = int(sys.argv[3])
    else:
        MIN_STREAK = 10

    sys.stderr.write( "Parsing...\n" )
    
    raw_data = [ [float(i) for i in line.split()]
                 for line in sys.stdin.readlines() ]
    
    # to rotate the plot
    data = [ [ j[i] for j in raw_data ]
             for i in range(len(raw_data[0])-1,-1,-1) ]

    ### Selecting only those points with value outside TOLERANCE

    sys.stderr.write( "Scanning...\n" )
    points = [ [0 for i in col] for col in data ]
    points_copy = [ [0 for i in col] for col in data ]
    linenum = 0
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

    ### Rerotate the image
    sys.stderr.write( "Printing...\n" )
    for i in range(len(points)):
        for j in range(len(points[0])):
            if points[i][j]==1 or points_copy[i][j]==1:
                print i,j
        print

if __name__=="__main__":
    main()
