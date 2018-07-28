mkdir dif
python sPlots.py
total=$(grep Nt constants.dat | cut -d "." -f 1 | cut -d " " -f 2 )
dens="python gif.py 0.14 "
acce="python gif.py 0.14 "
phase="python gif.py 0.14 "
pot="python gif.py 0.14 "
#echo $com
for i in $(seq 0 $(($total-1)))
do
#    echo $i
    phase+="dif/phase$i.png "
    dens+="dif/density$i.png "
    pot+="dif/potential$i.png "
    acce+="dif/acce$i.png "
done
eval $phase
eval $dens
eval $pot
eval $acce
