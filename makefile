all: grid.dat .PHONY

grid.dat: simulacion.c
	@echo Preparando y corriendo la simulación:
	mkdir datFiles
	mkdir images
	gcc -o heh.a -I/usr/local/include -L/usr/local/lib/ simulacion.c -lfftw3l -lm
	./heh.a
	python simPlots.py

.PHONY: clean
clean:
	@echo eliminando archivos temporales
	rm density.dat grid.dat heh.a oI.dat oR.dat outF0.dat outF1.dat inF.dat acce.dat
