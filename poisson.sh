make -f clean
mkdir datFiles
mkdir images
gcc -no-pie -o heh.a -I/usr/local/include -L/usr/local/lib/ poissonSolver.c -lfftw3 -lm 
./heh.a
rm heh.a
cd datFiles
cp density0.dat potential0.dat potential1.dat fdens0.dat fdens1.dat ~/2018/Varios
cd ..

