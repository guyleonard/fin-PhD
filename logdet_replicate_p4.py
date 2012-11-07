#!/usr/bin/python

# import all of p4 module
from p4 import *
import sys

# load alignment from argument when invoking
read(sys.argv[1])
in_alignment = var.alignments[0]

# output filename
output_all_trees_nexus=sys.argv[1]+"_all_trees_file.nxs"
consensus_tree_output=sys.argv[1]+"_consensus_tree.nxs"

# init lists
bootstrapped_alignments=[]
log_det_distance_matrices=[]
nj_trees=[]

# init taxnames
species_list=in_alignment.taxNames

# loop 1000 times
for index in range(1000):

    print "\n\nReplicate",index,":"

    # bootstrap replicate alignment
    bootstrapped_alignments.append(in_alignment.bootstrap())
    print "Bootstrap generated"

    # make logdet distance matrix for each replicate
    log_det_distance_matrices.append(bootstrapped_alignments[index].logDet(doPInvarOfConstants=True, missingCharacterStrategy='fudge', nonPositiveDetStrategy='invert'))
    print "LogDet distance matrix generated"

    # build nj tree for each distance matrix
    nj_trees.append(log_det_distance_matrices[index].bionj())
    print "NJ tree generated"

    # write tree to nexus for safe keeping
    nj_trees[index].writeNexus(fName=output_all_trees_nexus, append=1)


# make trees object containing all trees
all_trees_as_trees_obj=Trees(nj_trees,taxNames=species_list)

# make tree partition object from all trees
all_tree_partitions_as_tp_object=TreePartitions(all_trees_as_trees_obj)

# build consenus tree from tree partition object
consensus_tree = all_tree_partitions_as_tp_object.consensus()

# transfer support values to tree from support object
for n in consensus_tree.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)

# output consensus tree to stdout
print "\n\n\n\n\n\n Consensus Tree: \n\n\n\n"
consensus_tree.draw()

# print consensus tree to file
consensus_tree.writeNexus(consensus_tree_output)


