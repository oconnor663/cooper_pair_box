#! /usr/bin/python

import sys

def main():
    sys.stderr.write( "Getting data...\n" )
    raw_data = []
    for line in sys.stdin:
        raw_data.append([ float(j) for j in line.split() ])
        if len(raw_data)%100==0:
            sys.stderr.write( "Line %i\n" % len(raw_data) )
    # rotate without reversing
    sys.stderr.write( "Rotating...\n" )
    data = [ [ i[j] for i in raw_data ] for j in range(len(raw_data[0])-1,-1,-1) ]

    sys.stderr.write( "Outputting...\n" )
    for i in range(len(data)):
        if i%100==0:
            sys.stderr.write( "Column %i\n" % i )
        for j in range(len(data[0])):
            print "%i %i %f" % (i,j,data[i][j])
        print

main()
