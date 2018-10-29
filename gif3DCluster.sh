#!/bin/bash
total=$(grep Nt images/dats/datFiles/constants.dat | cut -d "." -f 1 | cut -d " " -f 2 )
echo $total
time="0.2 "
densXY="python gif.py "+$time
for i in $(seq 0 1 $(($total-1)))
do
#    echo $i
    densXY+="images/images/XY$i.png "
done
eval $densXY

