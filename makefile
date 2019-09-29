all: run plots gifs

clean:
	@echo Deleting output files from former execution.
	rm -rf ./datFiles
	rm -rf ./images
	rm -rf ./evolution

run: simulation.c clean
	@echo Creating output folders.
	mkdir datFiles
	mkdir images
	mkdir evolution
	@echo Compiling and executing.
	gcc -no-pie -o run.a -I/usr/local/include -I${HOME}/gitStuff/PoisFFT/src -L/usr/local/lib/ -L${HOME}/gitStuff/PoisFFT/lib/gcc  simulation.c -Wall -lpoisfft -lm  -lfftw3 -lfftw3f -lfftw3_omp -fopenmp
	touch run
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
	rm run

gifs: plots
	bash gif.sh
	rm plots

pyt: datFiles
	@echo drawing plots.
	python simPlots.py
	bash gif.sh

compare: 
	cp col/constants.dat constants.dat
	bash compare.sh
	rm constants.dat
