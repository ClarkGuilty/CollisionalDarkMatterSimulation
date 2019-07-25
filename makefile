all: clean data plots gifs .PHONY

clean:
	@echo Deleting output files from former execution.
	rm -rf ./datFiles
	rm -rf ./images


data: simulation.c
	@echo Creating output folders.
	mkdir datFiles
	mkdir images
	@echo Compiling and executing.
	gcc -no-pie -o run.a -I/usr/local/include -L/usr/local/lib/ simulation.c -Wall -lfftw3 -lm 
	./run.a
	rm run.a	

plots: data
	python simPlots.py

gifs: plots
	bash gif.sh

pyt: datFiles
	@echo drawing plots.
	python simPlots.py
	bash gif.sh

compare: 
	cp col/constants.dat constants.dat
	bash compare.sh
	rm constants.dat
