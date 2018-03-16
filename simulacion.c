/*
Javier Alejandro Acevedo Barroso
Primer Bosquejo. 1D con método de fourier.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Constantes de la simulación.
#define PI 3.14159265359

//Valores límites para la posición.
#define Xmin -1
#define Xmax 1

//Valores límites para la velocidad
#define Vmin -1
#define Vmax 1

#define Nx 2048
#define Nv 2048


//Arreglos
double phase[Nx][Nv];
double density[Nx];

//Variables
int i;
int j;
int k;

double dx = (Xmax-Xmin)*1.0/Nx;
double dv = (Vmax-Vmin)*1.0/Nv;

FILE *constantes;

//Métodos
void printPhase(char *name);
double gaussD(double x, double v, double sx, double sv, double amplitude);
double gaussD1(double x, double v, double sx, double sv);
double calDensity();
void printDensity(char *name);
void printConstant(char *name, double value);


int main()
{
	constantes = fopen("constants.dat","w+");
	printConstant("Xmin",Xmin);
	printConstant("Xmax",Xmax);
	printConstant("Vmin",Vmin);
	printConstant("Vmax",Vmax);
	printConstant("Nx",Nx);
	printConstant("Nv",Nv);
	double x;
	double v;
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nv;j+=1){
			x = Xmin*1.0+dx*i;
			v = Vmin*1.0+dv*j;
			phase[i][j] = gaussD(x,v,4,0.08, 0.08);
				}
			}
	printPhase("grid.dat");
	double mass = calDensity();
	printf("%f\n",mass);
	printDensity("density.dat");

	fclose(constantes);
	return 0;

}



void printPhase(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nv;j+=1){
			fprintf(output,"%f ", phase[i][j]);
				}
		fprintf(output,"\n");
			}
	fclose(output);

}

void printDensity(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",density[i]);
			}
	fclose(output);
}

void printConstant(char *name, double value)
{
    fprintf(constantes, "%s", name);
    fprintf(constantes, " %f\n", value);

}


double gaussD(double x, double v, double sx, double sv, double amplitude)
{
	double ex = -x*x/(2.0*sx*sx)-v*v/(2.0*sv*sv);
	return amplitude*exp(ex);

}

double gaussD1(double x, double v, double sx, double sv)
{
	double ex = -x*x/(2.0)-v*v/(2.0);
	double a = (1.0/(2.0*PI));
	return (1.0/Nx)*a*a*exp(ex);

}

double calDensity()
{
	double mass = 0;
	for(i = 0;i<Nx;i+=1){
		density[i]=0;
		for(j=0;j<Nv;j+=1){
				density[i] += phase[i][j]*dv;
			}
		mass += density[i];

		}
	return mass;
}


























