all: outputGo data plots gifs .PHONY

outputGo:
	@echo Deleting output files from former execution.
	rm -rf ./datFiles
	rm -rf ./images


data: simulacion.c
	@echo Creating output folders.
	mkdir datFiles
	mkdir images
	@echo Compiling and executing.
	gcc -no-pie -o heh.a -I/usr/local/include -L/usr/local/lib/ simulacion.c -Wall -lfftw3 -lm 
	./heh.a
	
plots: 
	python simPlots.py

gifs: images
	bash gif.sh

.PHONY: clean

clean:
	@echo Removing temporal files
	rm heh.a

pyt: 
	@echo drawing plots.
	python simPlots.py
	bash gif.sh
