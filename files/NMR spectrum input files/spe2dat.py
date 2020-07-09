#!/usr/bin/env python

import sys

print "Convert Simpson SPE files to plottable format"
spe = raw_input("SPE file: ")
larmor = float(raw_input("Larmor freq. (MHz): "))

try:
    spef = open(spe)
except:
    print("cannot open file", spe)
    sys.exit(1)

# parse SPE header
err = 0
if spef.readline()[0:4] != "SIMP": err = 1
np = int(spef.readline()[3:])
sw = float(spef.readline()[3:])
if spef.readline()[0:8] != "TYPE=SPE": err = 1
if spef.readline()[0:4] != "DATA": err = 1

if err == 1:
    print("file format error")
    sys.exit(1)

# read SPE data
data = []
for i in xrange(np):
    a, b = spef.readline().split()
    data.append(float(a))
spef.close()


# write data
out = open(spe+".dat", "wt")
for i in xrange(np):
    freq = i*sw/np - sw/2.0
    print >>out, freq/larmor, data[i]
out.close()

print("done")


