#! /usr/bin/python

import sys

def main():
    input = sys.stdin.readlines()
    raw_data = [ [ float(j) for j in i.split() ] for i in input ]
    # rotate without reversing
    data = [ [ i[j] for i in raw_data ] for j in range(len(raw_data[0])-1,-1,-1) ]

    for i in range(len(data)):
        for j in range(len(data[0])):
            print "%i %i %f" % (i,j,data[i][j])
        print

main()
