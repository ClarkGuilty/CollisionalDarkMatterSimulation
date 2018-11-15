mkdir dif
echo "Generando Plots"
python plotsCompare1D.py
echo "Generando gif"
total=$(grep Nt constants.dat | cut -d "." -f 1 | cut -d " " -f 2 )
dens="python gif.py 0.4 "
#acce="python gif.py 0.5 "
#pot="python gif.py 0.5 "
phase="python gif.py 0.4 "
#echo $com
for i in $(seq 0 1)
do
#    echo $i
    phase+="dif/phase$i.png "
    dens+="dif/density$i.png "
#    pot+="dif/potential$i.png "
#    acce+="dif/acce$i.png "
done
for i in $(seq 2 2 $(($total-2)))
do
#    echo $i
    phase+="dif/phase$i.png "
    dens+="dif/density$i.png "
#    pot+="dif/potential$i.png "
#    acce+="dif/acce$i.png "
done
eval $phase
eval $dens
#eval $pot
#eval $acce
