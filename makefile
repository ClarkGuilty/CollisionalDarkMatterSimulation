all: grid.dat .PHONY

grid.dat: simulacion.c
	@echo corriendo la simulación
	gcc -o heh.a -I/usr/include -L/usr/lib/x86_64-linux-gnu simulacion.c -lfftw3 -lm
	./heh.a
	python simPlots.py

.PHONY: clean
clean:
	@echo eliminando archivos temporales
	rm density.dat grid.dat heh.a oI.dat oR.dat outF0.dat outF1.dat inF.dat acce.dat grid1.dat grid2.dat
