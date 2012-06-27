#!/bin/bash
if [$1 = ]; then
    echo "no folder name"
    exit 0  
fi

folder=$1

mkdir $folder

cd $folder

mkdir 1

cd 1

mkdir 1.sequence 2.refgen 3.alignmnet 4.mask 5.model 6.tree 7.treenamer

cd 6.tree 

mkdir phyml mrbayes raxml

exit 0
