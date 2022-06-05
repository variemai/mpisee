#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    print("Input mpiP profile output")
    exit(0)

with open(sys.argv[1],'r') as f:
    lines = f.readlines()
    found = 0
    counter = 1
    for line in lines:
        if line.find("Callsite Time") != -1:
            found = 1
        if line.find("Callsite Message") != -1:
            break
        if ( found and line.find("Bcast") != -1 ):
            counter += 1
    print ( counter )
