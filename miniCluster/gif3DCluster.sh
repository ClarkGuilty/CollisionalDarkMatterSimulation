#!/bin/bash
total=$(grep Nt datFiles/constants.dat | cut -d "." -f 1 | cut -d " " -f 2 )
#total=50
echo $total
time="0.2 "
densXY="python gif.py "+$time
for i in $(seq 0 2 $(($total-1)))
do
#    echo $i
    densXY+="images/XY$i.png "
done
eval $densXY

