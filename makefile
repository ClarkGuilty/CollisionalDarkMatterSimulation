all: grid.dat .PHONY

grid.dat: simulacion.c
	@echo corriendo la simulación
	gcc -o heh.a simulacion.c -lm
	./heh.a
	python simPlots.py

.PHONY: clean
clean:
	@echo eliminando archivos temporales
	rm density.dat grid.dat heh.a
