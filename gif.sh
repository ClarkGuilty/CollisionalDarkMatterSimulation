#!/bin/bash
total=$(grep Nt constants.dat | cut -d "." -f 1 | cut -d " " -f 2 )
com="python gif.py 0.07 "
#echo $com
for i in $(seq 0 $(($total-1)))
do
#    echo $i
    com+="images/phase$i.png "
done
eval $com
