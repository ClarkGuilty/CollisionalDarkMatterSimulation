all: clean run.a  plots gifs

clean:
	@echo Deleting output files from former execution.
	rm -rf ./datFiles
	rm -rf ./images
	rm -rf ./evolution

run.a: simulation.c
	@echo Creating output folders.
	mkdir datFiles
	mkdir images
	mkdir evolution
	@echo Compiling and executing.
	gcc -no-pie -o run.a -I/usr/local/include -L/usr/local/lib/ simulation.c -Wall -lfftw3 -lm 
	./run.a
	rm run.a	

tests: datFiles/constants.dat testPlots.py
	cp fourierEvolution.dat ~/Paper-CollDM/plots
	cp datFiles/constants.dat ~/Paper-CollDM/plots
	python ~/Paper-CollDM/plots/plots.py
	python testPlots.py
	cp images/powerSeries8.png ~/Paper-CollDM/plots
	cp images/powerSeries16.png ~/Paper-CollDM/plots
	cp images/powerSeries24.png ~/Paper-CollDM/plots
	cp images/powerSeries32.png ~/Paper-CollDM/plots


plots: datFiles/constants.dat simPlots.py
	python simPlots.py

gifs: plots
	bash gif.sh
	rm plots

poissonTest: clean
	@echo Creating output folders.
	mkdir datFiles
	mkdir images
	mkdir evolution
	@echo Compiling and executing.
	gcc -no-pie -o run.a -I/usr/local/include -I/home/hiparco/gitStuff/PoisFFT/src -L/usr/local/lib/ -L/home/hiparco/gitStuff/PoisFFT/lib/gcc  simulationPoiss.c -Wall -lpoisfft -lm  -lfftw3 -lfftw3f -lfftw3_omp -fopenmp
	./run.a
	rm run.a



pyt: datFiles
	@echo drawing plots.
	python simPlots.py
	bash gif.sh

compare: 
	cp col/constants.dat constants.dat
	bash compare.sh
	rm constants.dat
