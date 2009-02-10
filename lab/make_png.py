#! /usr/bin/python

import os,sys

def main():
    out = os.popen("gnuplot","w")
    out.write( "set terminal png size 300,500\n"
               "set output 'out.png'\n"
               "plot '-' with points lt 3 pt 0\n" )
    line = sys.stdin.readline()
    while line:
        out.write(line)
        line = sys.stdin.readline()
    out.write("e\n")
    out.close()

main()
