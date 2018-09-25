#!/bin/bash
total=$(grep Nt constants.dat | cut -d "." -f 1 | cut -d " " -f 2 )
dens="python gif.py 0.2 "
acce="python gif.py 0.2 "
phasex="python gif.py 0.2 "
phasey="python gif.py 0.2 "
pot="python gif.py 0.2 "
for i in $(seq 0 2 $(($total-1)))
do
#    echo $i
    phasex+="images/gridx$i.png "
    phasey+="images/gridy$i.png "
    dens+="images/density$i.png "
#    pot+="images/potential$i.png "
#    acce+="images/acce$i.png "
done
eval $phasex
eval $phasey
eval $dens
#eval $pot
#eval $acce
