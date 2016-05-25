#!/bin/bash

Tempurature=(1 1.5 2 2.5 3 3.5 4 4.5)
# Input Tempurature you want to test, use space to separate different value
for a in ${Tempurature[@]}; do
    echo Tempurature: $a    Outputfile: isingT$a.txt
    ./ising.o $a  >   isingT$a.txt
done
