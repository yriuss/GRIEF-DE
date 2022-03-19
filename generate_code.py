#! /usr/bin/env python

import numpy
import sys

tests = numpy.loadtxt(sys.argv[1], dtype=int)
if len(sys.argv) > 2:
    bytes = int(sys.argv[2])
else:
    bytes = tests.shape[0] / 8

print ('// Code generated with \'$ %s\'' % ' '.join(sys.argv))
print ('#define SMOOTHED(y,x) smoothedSum(sum, pt, y, x)')
for byte in range(bytes):
    print ('    desc[%d] =' % byte,)
    for i in range(8):
        x1, y1, x2, y2 = tests[byte*8 + i, :]
        print ("((SMOOTHED(individual[" + str(byte*8 + i) + "][" + str(0) + "], individual[" + str(byte*8 + i) + "][" + str(1) + "]) < SMOOTHED(individual[" + str(byte*8 + i) + "][" + str(2) + "], individual[" + str(byte*8 + i) + "][" + str(3) + "])) << " + str(7-i) + ")")
        if i == 7:
            print (';')
        else:
            print ('+',)
print ('#undef SMOOTHED')
