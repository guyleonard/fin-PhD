#!/usr/bin/python
# Author:Finlay Maguire
# Date:09/11/2012

# Usage: python logdet_replicate_p4.py <INPUT ALIGNMENT> <DEBUG STATUS 0/1>

# import all of p4 module
from p4 import *
import sys

# set debug status
debug=int(sys.argv[2])

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

    if debug==1:
        print "\n\n-----ALIGNMENT-----\n\n"
        print bootstrapped_alignments[index]
        bootstrapped_alignments[index].writeNexus(fName=None)

    print "Bootstrap generated"

    # make logdet distance matrix for each replicate
    ld_mat=bootstrapped_alignments[index].logDet(correction='TK02', doPInvarOfConstants=True, missingCharacterStrategy='fudge', nonPositiveDetStrategy='invert')
    log_det_distance_matrices.append(ld_mat)

    if debug==1:
        print "\n\n-----MATRIX-----\n\n"
        print log_det_distance_matrices[index]
        log_det_distance_matrices[index].writeNexus(fName=None)

    print "LogDet distance matrix generated"

    # build nj tree for each distance matrix
    tree=log_det_distance_matrices[index].bionj()
    nj_trees.append(tree)

    if debug==1:
        print "\n\n-----PHYLOGENY----\n\n"
        print nj_trees[index]
        nj_trees[index].writeNexus(fName=None)

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

