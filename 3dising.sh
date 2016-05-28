#!/bin/bash

Tempurature=(1.0 2.0 3.0 4.0 5.0 6.0)
LATTICE_LENGTH=(20)
# Input Tempurature you want to test, use space to separate different value

for l in ${LATTICE_LENGTH[@]}; do
    for t in ${Tempurature[@]}; do
        echo Tempurature: $t lattice_length $l   Outputfile: ising3dT${t}L${l}.txt
        ./src/ising3d.o $t  $l >   data/3d/ising3dT${t}L${l}.txt
    done
done
