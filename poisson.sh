make -f clean
mkdir datFiles
mkdir images
gcc -no-pie -o heh.a -I/usr/local/include -L/usr/local/lib/ poissonSolver.c -Wall -lfftw3 -lm 
./heh.a
rm heh.a

