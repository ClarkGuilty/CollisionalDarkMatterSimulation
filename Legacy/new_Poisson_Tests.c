/*
Written by Javier Alejandro Acevedo Barroso
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>

//Constants of the simulation.
#define PI 3.14159265359

//Extreme values for velocity and position.
#define Xmin -1.0
#define Xmax 1.0
#define Vmin -1.0
#define Vmax 1.0
#define scale 1 //scale in order to get better graphics. Better left at 1 after all.

//Grid size.
#define Nx 16584
#define Nv 2048

//Int constants for comparision.
#define aKpc 18
#define aSegundos 14
#define aByear 4
#define aMasasSol 5
#define aMetros 464
 
#define GAUSS -127
#define JEANS -137
#define BULLET -147

//Relaxation time. 0 means solving the Vlasov equation.
//#define TAU 8972
#define TAU 0

//Units of the simulation. This particular set corresponds to a galactic scale.
#define mParsecs 35e-3  //How many mpc are equivalent to one spatial unit.
#define solarMases 1e12 //How many solar masses are equivalent to one mass unit.
#define fracT0 4e-3     //What fraction of the age of the universe is equivalent to one time unit.
//#define G 0.959572 //Gravitational constant in this units. It is calculated with sPlots.py
#define G 0.031830 

//Set of units for a galaxy cluster scale.
//#define mParsecs 5
//#define solarMases 1e15
//#define fracT0 2e-1
//#define G 0.2729448134597113


//Arrays
double phase[Nx][Nv] = {0};

double phaseOld[Nx][Nv] = {0};
double phaseTemp[Nx][Nv] = {0};
double *energy;
double *velocity;
double *density;
double *pot;
double *acce;
int initCon;
double totalMass;

//Variables and parameters
int i;
int j;
int k;
int l;
int i2;
int j2;

//Size of the grid
double Lx = Xmax - Xmin;
double Lv = Vmax - Vmin;
double dx = (Xmax-Xmin)*1.0/Nx;
double dv = (Vmax-Vmin)*1.0/Nv;

//Size of a timestep and number of timesteps.
double dt = 0.2; 
int Nt = 100;

//File with parameters of the simulation.
FILE *constants;


//Methods
double testPot(double x, double sx, double amplitude);
double testDens(double x, double sx, double amplitude);
void printPhase(char *name);
double gaussD(double x, double v, double sx, double sv, double amplitude);
double calDensity();
void printDensity(char *name);
void printConstant(char *name, double value);
double giveDensity(int l);
void potential();
double calcK2(int j2);
double calcK2T(int j2);
double convert(double valor, int unidad);
void calAcceG();
void printAcce(char *name);
double newij(int iin, int jin);
void step();
int mod(int p, int q);
void printPot(char *name);
double einasto(double x, double v, double sx, double sv, double amplitude);
double jeans(double x, double v, double rho, double sigma, double A, double k);
double feq(int ipos, int jvel);
double feq2(int ipos, int jvel);
double givePos(int ito);
double giveVel(int jto);
double collision(int icol, int jcol, double tau);
void collisionStep();
double newijCol(int iin, int jin);
double bulletC(double x, double v, double sx1, double sx2, double sv, double amplitude1,double amplitude2);

int main()
{
    //dt = dt*dx/dv;

    //Initialization of arrays.
    energy = malloc((sizeof(double)*Nx));
    velocity = malloc((sizeof(double)*Nx));
    density = malloc((sizeof(double)*Nx));
    acce = malloc((sizeof(double)*Nx));
    pot = malloc((sizeof(double)*Nx));
    
	double x;
    
    //Choosing what initial conditions to simulate.
    //initCon = GAUSS;
    initCon = JEANS;
    //initCon = BULLET;
    

    
    //Initialization parameters
    //Poisson solver test.
    double sx = 0.1;
    double ampl = 40.0;
    
    
    //Poisson test.
    
	for(i=0;i<Nx;i+=1) {
                x = Xmin*1.0+dx*i;
                density[i] = testDens(x,sx,ampl);
                pot[i] = testPot(x,sx,ampl);
        }

		
	
    
	printPot("./testPotential.dat");
	printDensity("./testDensity.dat");
    
    
    //Solves Poisson equation.
    potential();
    
    printPot("./solvedPotential.dat");
    
    
	
	return 0;

}



//gaussian potential.
double testPot(double x, double sx, double amplitude)
{
    return amplitude*exp(-x*x/(sx*sx));
    
}

//Corresponding density for a gaussian potential.
double testDens(double x, double sx, double amplitude)
{
    return 2*(2*x*x-sx*sx)*amplitude*exp(-x*x/(sx*sx))/(sx*sx*sx*sx);
}

//Exports the phase space grid to text. char name is the name of the output file.
void printPhase(char *name)
{
	FILE *output = fopen(name, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Nv+1;j+=1){ 
          //      printf("ignorarPrimero\n");
			fprintf(output,"%f ", convert(phase[i][Nv-j], aMasasSol)/convert(1.0,aKpc)/(convert(1.0,aKpc)*3.0857e+19)* convert(1.0,aSegundos)); //Imprime en Masas solares /kpc / (km/s)
            //fprintf(output,"%f ",phase[i][Nv-j]);
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);

}
//Exports the density array to text. char name is the name of the output file.
void printDensity(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",convert(density[i],aMasasSol)/convert(1,aKpc)); //Imprime en Masas solares / kiloparsec.
			}
	fclose(output);
}

//Calcula la densidad, la velocidad macroscópica (u) y la energía libre (e) y los carga en los correspondientes arreglos.
//Integrates the phase space with regards to velocity in order to obtain density, average velocity (u) and local free energy (e).
double calDensity()
{
	double mass = 0;
	for(i = 0;i<Nx;i+=1){
		density[i]=0;
        velocity[i] = 0;
        energy[i] = 0;
		for(j=0;j<Nv;j+=1){
				density[i] += phase[i][j]*dv;
                velocity[i] += phase[i][j]*dv*giveVel(j);
        }
        for(j=0;j<Nv;j+=1){
            energy[i] += phase[i][j]*dv*pow(giveVel(j) - velocity[i], 2)/2.0;
        }
        velocity[i] = velocity[i] / density[i];
        energy[i] = energy[i] / density[i];
		mass += density[i]*dx;
		}
	return mass;
}

//Solves the Poisson equation using Fourier Method. Updates the potential array.
void potential()
{

    fftw_complex *in, *out, *inR, *mem;
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    inR=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    mem=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);


    fftw_plan pIda;
    pIda = fftw_plan_dft_1d(Nx, in, out,FFTW_FORWARD, FFTW_MEASURE);

    //Cargar densidad en in y borra out:
    for(i=0;i<Nx;i+=1){

        in[i] = giveDensity(i) - totalMass/(Xmax-Xmin);
        inR[i] = -1.0;
        out[i] = 0;
    }
    fftw_execute(pIda);

    //Guarda out en mem.
    for(i=0;i<Nx;i+=1){
        mem[i] = out[i];
    }

    pIda = fftw_plan_dft_1d(Nx, out, inR, FFTW_BACKWARD, FFTW_MEASURE); //Se debe usar el mismo plan sí o sí al parecer.
    
    //Devuelve carga a out Î(Chi).
    out[0] = -mem[0];
    for(i=1;i<Nx;i+=1){
      out[i] = -mem[i]*calcK2(i);
    //out[i] = mem[i]; //Descomentar esta línea para obtener la distribucion original.
    }
    fftw_execute(pIda);

    
    for(i=0;i<Nx;i+=1){
        pot[i] = creal(inR[i]/Nx);
    }
}


double calcK2(int j2)
{
    double rta;
    if( ( (j2 == 0)   )  ){
        return 0;
    }
    if(j2<Nx/2){
        rta = 2*PI*j2/(Xmax-Xmin);
    }
    if(j2>=Nx/2){
        rta = -2*PI*(Nx-j2)/(Xmax-Xmin);
    }    
    return pow(1.0/rta,2);
}

double calcK2_old(int j2)
{
    double rta;
    if( ( (j2 == 0)   )  ){
        return 0;
    }
    if ( (j2 == Nx/2+1) ) {
        return 0;
    }
    
    if(j2<Nx/2+1){
        rta = PI*j2;
    }
    if(j2>=Nx/2+1){
        rta = -PI*(Nx-j2);
    }    
    return pow(1.0/rta,2);

}


double calcK2_backup(int j2)
{
    double rta;
    if( ( (j2 == 0)   )  ){
        return 0;
    }
    if ( (j2 == Nx/2+1) ) {
        return 0;
    }
    
    if(j2<Nx/2+1){
        rta = 2*PI*j2/(Xmax-Xmin);
    }
    if(j2>=Nx/2+1){
        rta = -2*PI*(Nx-j2)/(Xmax-Xmin);
    }    
    return pow(1.0/rta,2);
}

double calcK2T(int j2)
{
    double rta;
    if( ( (j2 == 0))){
        return 0;
    }
//    if ( (j2 == Nx/2+1) ) {
        //return 0;
    //}    
    if(j2<Nx/2){
        rta = 2*sin(PI*j2/(Nx))/dx;
    }
    if(j2>=Nx/2){
        //rta = -2*PI*(Nx-j2)/(Xmax-Xmin);
        rta = 2*sin((Nx-j2)*PI/(Nx))/dx;
        //rta = 2*sin(PI*(Nx-j2)/(Nx-1))/(Xmax-Xmin);
    }    
    return pow(1.0/rta,2);
}

//Método para evitar efectos misticos en la memoria.
double giveDensity(int l)
{
    double rta = density[l];
    return rta;
}

//Convierte unidades de la simulación a masas solares, metros, o segundos.
double convert(double valor, int unidad )
{
    double conx0 = 3.0857e+22; //un megaparsec en metros
    double cont0 = 13.772*1000000000; //Edad del universo
    double cont0s = cont0*365.24*24*60*60;

    if(unidad == aMasasSol){
        return valor * solarMases;
    }
    if(unidad == aKpc){
        return valor * mParsecs * 1000.0; //1000 = 1megaparsec en kiloparsec
    }
    if( unidad == aByear){
        return valor*13.772*fracT0;
    }
    if( unidad == aSegundos){
        return valor*cont0s*fracT0;
    }
    if(unidad == aMetros){
        return valor * mParsecs * conx0;
    }
    return -1;
}

//Deriva el potential y carga la aceleración gravitacional en el arreglo acce.
void calAcceG()
{
    for(i = 0; i<Nx ; i +=1){
    acce[i] =  (-pot[mod(i-2,Nx)] + 8*pot[mod(i-1,Nx)]-8*pot[mod(i+1,Nx)]+pot[mod(i+2,Nx)])/(12*dx);
    }
}

//Imprime el arreglo Acce.
void printAcce(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",convert(acce[i], aKpc)/pow(convert(1.0, aByear)*1000,2)); //Imprime en kpc / (mAños)^2;
			}
	fclose(output);
}

//Imprime el arreglo Pot
void printPot(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",pow(convert(pot[i], aMetros)/convert(1.0, aSegundos),2)/pot[i]); //Imprime potential en J/kg
			}
	fclose(output);
}

//Calcula la nueva posición de un elemento de la grilla tras el streaming. Version que discretiza los cambios y no el nuevo valor.
double newij(int iin, int jin)
{
        double x = Xmin*1.0+dx*iin; //Inicialización

        double v = acce[iin]*dt;
        double dj = v/dv;
        dj = (int)dj;
        j2 = jin+dj;

        if(j2 < 0 || j2 >= Nv) return -1;
        v = Vmin*1.0+dv*j2;
        x = v*dt;
        double di = x/dx*scale;
        di = (int) di;

        i2 = iin + di;
        i2 = mod(i2,Nx);
    return 0;
}

//Calcula un paso. Guarda una copia del phase actual en phaseOld. Actualiza phase. k,i son x. j,l son v.
void step()
{
	for(k = 0; k<Nx; k++){
		for(l= 0; l<Nv; l++){
			if(newij(k,l) ==0){
				phaseTemp[i2][j2] += phase[k][l];
			}
		}
	}

	for(i = 0; i<Nx; i++){
		for(j= 0; j<Nv; j++){
			phase[i][j] = phaseTemp[i][j];
			phaseTemp[i][j] = 0;
		}
	}
}

void collisionStep()
{
    	for(k = 0; k<Nx; k++){
		for(l= 0; l<Nv; l++){
			if(newijCol(k,l) ==0){
                phaseTemp[i2][l] += collision(k,l,TAU) + phase[k][l] ;//+ dt*feq2(k,l)*acce[k]*(giveVel(l)-velocity[k])/energy[k];
			}
		}
	}

	for(i = 0; i<Nx; i++){
		for(j= 0; j<Nv; j++){
			phase[i][j] = phaseTemp[i][j];
			phaseTemp[i][j] = 0;
		}
	}
}

//Calcula el cambio en x. Ignora j.
double newijCol(int iin, int jin)
{
        double x = Xmin*1.0+dx*iin; //Inicialización


        double v = acce[iin]*dt;
        j2 = jin; 

        if(j2 < 0 || j2 >= Nv) return -1;
        v = giveVel(j2);
        
        x = v*dt*scale;
        double di = x/dx;
        di = (int) di;

        i2 = iin + di;
        i2 = mod(i2,Nx);
    return 0;
}

//Calcula la contribución colisional en phase[icol][jcol] con un Tau dado.
double collision(int icol, int jcol, double tau)
{
    if(TAU==0) return 0;
    double df = (feq(icol,jcol) - phase[icol][jcol])/(1.0*tau);    
    return df;
}

double feq(int ipos, int jvel)
{
    double ex = -1.0*pow(giveVel(jvel)-velocity[ipos],2)/(2.0*energy[ipos]);
    double other = density[ipos] / sqrt(2.0*PI*energy[ipos]);
    return other * exp(ex);    
}

double feq2(int ipos, int jvel)
{
    double ex = -1.0*pow(giveVel(jvel),2)/(2.0*energy[ipos]);
    double other = density[ipos] / sqrt(2*PI*energy[ipos]);
    double lowMach = 1.0 + giveVel(jvel)*velocity[ipos]/energy[ipos] + pow(giveVel(jvel)*velocity[ipos],2)/(2.0*energy[ipos]) - pow(velocity[ipos],2)/(2.0*energy[ipos]);
    return other * exp(ex)* lowMach;
}


double givePos(int ito)
{
    return Xmin*1.0+dx*ito;
}

double giveVel(int jto)
{
    return Vmin*1.0+dv*jto;
}

//Observando que el m = p % q es negativo si p<0 y q>0, se define una función de módulo con rango de 0 a q-1.
int mod(int p, int q)
{
	p = p%q;
	if(p<0){
		 return p+q;
	}
	return p;

}

















