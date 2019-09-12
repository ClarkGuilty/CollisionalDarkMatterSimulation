all: clean run.a plots gifs

clean:
	@echo Deleting output files from former execution.
	rm -rf ./datFiles
	rm -rf ./images


run.a: simulation.c
	@echo Creating output folders.
	mkdir datFiles
	mkdir images
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

pyt: datFiles
	@echo drawing plots.
	python simPlots.py
	bash gif.sh

compare: 
	cp col/constants.dat constants.dat
	bash compare.sh
	rm constants.dat
