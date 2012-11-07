#!/usr/bin/python

from p4 import *
import sys

# load alignment from argument when invoking
read(sys.argv[1])
in_alignment = var.alignments[0]

# initialise lists
bootstrapped_alignments=[]
log_det_distance_matrices=[]
nj_trees=[]


#loop 1000 times
for index in range(2):

    #bootstrap replicate alignment
    bootstrapped_alignments.append(in_alignment.bootstrap())
    print "Bootstrap",index,"generated"

    #make logdet distance matrix for each replicate
    log_det_distance_matrices.append(bootstrapped_alignments[index].logDet(doPInvarOfConstants=True, missingCharacterStrategy='fudge', nonPositiveDetStrategy='invert'))
    print "LogDet distance matrix",index,"generated"

    #nj tree each distance matrix
    nj_trees.append(log_det_distance_matrices[index].bionj)
    print "NJ tree",index,"generated"

print nj_trees

consensus_trees = Trees(nj_trees)

consensus_trees.writeNewick(fname="out")

