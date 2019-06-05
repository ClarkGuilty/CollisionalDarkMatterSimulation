#!/bin/bash
total=$(grep Nt datFiles/constants.dat | cut -d "." -f 1 | cut -d " " -f 2 )
dens="python gif.py 0.2 "
acce="python gif.py 0.2 "
phase="python gif.py 0.2 "
pot="python gif.py 0.2 "
for i in $(seq 0 2 $(($total-1)))
do
#    echo $i
    phase+="images/phase$i.png "
    dens+="images/density$i.png "
    pot+="images/potential$i.png "
    acce+="images/acce$i.png "
done
eval $phase
eval $dens
#eval $pot
#eval $acce
