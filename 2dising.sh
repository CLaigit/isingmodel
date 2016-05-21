#!/bin/bash

Tempurature=(1 1.5 2 2.5 3 3.5 4 4.5)
for a in ${Tempurature[@]}; do
    echo Tempurature: $a    Outputfile: isingT$a.dat
    ./ising $a  >   isingT$a.dat
done
