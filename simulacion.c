/*
Javier Alejandro Acevedo Barroso
Primer Bosquejo. 1D con método de fourier.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

//Constantes de la simulación.
#define PI 3.14159265359

//Valores límites para la posición y velocidad.
#define Xmin -1.0
#define Xmax 1.0
#define Vmin -0.5
#define Vmax 0.5
#define scale 0.5 //Re escala la velocidad para tener mejores gráficas.

//Tamaño del espacio.
#define Nx 2048
#define Nv 2048

//Constantes de unidades.
#define aKpc 18
#define aSegundos 14
#define aByear 4
#define aMasasSol 5
#define aMetros 464
 
#define GAUSS -127
#define JEANS -137
#define TAU 0


//Primer intento Via Láctea.
#define mParsecs 20e-3  //Cuántos megaparsecs equivalen a una unidad espacial.
#define solarMases 1e12 //Cuántas masas solares equivalen a una unidad de masa.
#define fracT0 3e-3     //Qué fracción de la edad del universo equivale a una unidad de tiempo
#define G 0.959572 //G en estas unidades. Se calcula con sPlots.py

//Unidades funcionales para clusters galácticos.
//#define mParsecs 5
//#define solarMases 1e15
//#define fracT0 2e-1
//#define G 0.2729448134597113


//Arreglos
double phase[Nx][Nv] = {0};

double phaseOld[Nx][Nv] = {0};
double phaseTemp[Nx][Nv] = {0};
double *energy;
double *velocity;
double *density;
double *pot;
double *acce;
int initCon;

//Variables
int i;
int j;
int k;
int l;
int i2;
int j2;
double Lx = Xmax- Xmin;
double Lv = Vmax- Vmin;
double dx = (Xmax-Xmin)*1.0/Nx;
double dv = (Vmax-Vmin)*1.0/Nv;

double dt = 0.5; 
int Nt = 100;
FILE *constantes;
void printPhase(char *name);
double gaussD(double x, double v, double sx, double sv, double amplitude);
double calDensity();
void printDensity(char *name);
void printConstant(char *name, double value);
double giveDensity(int l);
void potencial();
double calcK2(double j2);
double convertir(double valor, int unidad);
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

int main()
{
    //dt = dt*dx/dv;
    energy = malloc((sizeof(double)*Nx));
    velocity = malloc((sizeof(double)*Nx));
    density = malloc((sizeof(double)*Nx));
    acce = malloc((sizeof(double)*Nx));
    pot = malloc((sizeof(double)*Nx));
    

	constantes = fopen("constants.dat","w+");
	printConstant("Xmin",Xmin);
	printConstant("Xmax",Xmax);
	printConstant("Vmin",Vmin);
	printConstant("Vmax",Vmax);
	printConstant("Nx",Nx);
	printConstant("Nv",Nv);
	printConstant("Nt", Nt);
    
	double x;
	double v;
    
    //Variable que elige condición a simular.
    initCon = GAUSS;
    
    printConstant("InitCon", initCon);
    printConstant("TAU", TAU);
    
    //Gauss
    double vSx = 0.06;
    double vSv = 0.12;//Se reescala la velocidad al momento del drift, esto para mayor nitidez.
    double ampl = 1;
    
    //Jeans
    double rho = 0.1;
    double sigma = 0.1;
    double A = 0.9999;
    double k = 2* PI;
    
    
	for(i=0;i<Nx;i+=1) {
                x = Xmin*1.0+dx*i;
                //printf("x es %f en %d\n", x,i);
                //x = fabs(x);
                    for(j=0;j<Nv;j+=1){
                        v = Vmin*1.0+dv*j;
                        if(initCon == GAUSS)
                        {
                            phase[i][j] = gaussD(x,v,vSx,vSv,ampl);
                        }
                        if(initCon == JEANS)
                        {
                            phase[i][j] = jeans(x, v, rho, sigma, A, k);
                            //if(phase[i][j]<0) printf("(%d,%.2f) (%d,%.2f) = %.2f\n",i,j, givePos(i),giveVel(i), phase[i][j]);
                        }
                        phaseOld[i][j] = 0;
                        phaseTemp[i][j] = 0;
                    }
			}

		
	//printPhase("grid1.dat");
	double mass = calDensity();

	printf("Se simuló %f millones de años con %d pasos de %f millones de años cada uno\n", convertir(Nt*dt,aByear)*1000,Nt, convertir(dt,aByear)*1000);
	printf("La masa fue %f masas de la vía láctea\n",convertir(mass, aMasasSol)/1e12);
        printPhase("./datFiles/grid0.dat");
	printDensity("./datFiles/density0.dat");

	FILE *simInfo = fopen("./images/simInfo.dat","w+");
	printf("heh\n");
	fprintf(simInfo,"Para la simulación se utilizó las siguientes condiciones:\n");
	fprintf(simInfo,"x va de %.2f a %.2f , v va de %.2f a %.2f\n", Xmin,Xmax,Vmin,Vmax);
	if(initCon == GAUSS)
    {
        fprintf(simInfo,"Una distribución gaussiana centrada en 0 para el espacio de fase con (sx sv A)=\n");
        fprintf(simInfo,"(%.3f %.3f %.3f)\n", vSx, vSv, ampl);        
    }
    if(initCon == JEANS)
    {
        fprintf(simInfo,"Se simuló la inestabilidad de Jeans con (rho sigma A k)=\n");
        fprintf(simInfo,"(%.3f %.3f %.3f %.3f)\n", rho, sigma, A, k);
    }

	fprintf(simInfo,"Se simuló %d instantes con dt = %.3f\n", Nt,dt);
    
    potencial();
    printPot("./datFiles/potential0.dat");
    calAcceG();
    
    printAcce("./datFiles/acce0.dat");
    printf("G es %lf\n", G*1.0);


	for(int suprai = 1; suprai<Nt;suprai+=1){
        char *grid = (char*) malloc(200* sizeof(char));
        //printf("Error Mesage00\n");
		
		step();
        
        //Descomentar para versión colisional --
        if(TAU != 0){
        calDensity();  //La impresión calcula la densidad.
        collisionStep();
        }
        //--
		printf("%d %f\n",suprai,calDensity()*100/mass);
        
		sprintf(grid, "./datFiles/density%d.dat", suprai);
		printDensity(grid);

		potencial();
		sprintf(grid, "./datFiles/potential%d.dat", suprai);
		printPot(grid);
		calAcceG();
		sprintf(grid, "./datFiles/acce%d.dat", suprai);
		printAcce(grid);
                sprintf(grid, "./datFiles/grid%d.dat", suprai);
                printPhase(grid);
                free(grid);
                
	}


	fclose(constantes);
	fclose(simInfo);
	return 0;

}



//Imprime el espacio de fase con el String name como nombre.
void printPhase(char *name)
{
	FILE *output = fopen(name, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Nv+1;j+=1){ 
          //      printf("ignorarPrimero\n");
			fprintf(output,"%f ", convertir(phase[i][Nv-j], aMasasSol)/convertir(1.0,aKpc)/(convertir(1.0,aKpc)*3.0857e+19)* convertir(1.0,aSegundos)); //Imprime en Masas solares /kpc / (km/s)
            //fprintf(output,"%f ",phase[i][Nv-j]);
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);

}
//Imprime el arreglo density con el String name como nombre.
void printDensity(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",convertir(density[i],aMasasSol)/convertir(1,aKpc)); //Imprime en Masas solares / kiloparsec.
			}
	fclose(output);
}

//Imprime las constantes de la simulación para referencia fuera del C.
void printConstant(char *name, double value)
{
    fprintf(constantes, "%s", name);
    fprintf(constantes, " %f\n", value);

}

//Retorna el valor de la gaussiana para un x,v, sigma x, sigma v, y una amplitud dada.
double gaussD(double x, double v, double sx, double sv, double amplitude)
{
	double ex = -x*x/(2.0*sx*sx)-v*v/(2.0*sv*sv);
//	double ex = -x*x/(sx*sx)-v*v/(sv*sv);

	//return amplitude*exp(ex)/(2*PI*sx*sv);
	return amplitude*exp(ex);

}

//Retorna el valor correspondiente para la inestabilidad de Jeans. Con una densidad, una dispersión, una amplitud(0< A <=1) y un k dados.
double jeans(double x, double v, double rho, double sigma, double A, double k)
{
 return rho*pow(2*PI*sigma*sigma,-0.5)*exp(-v*v/(2*sigma*sigma))*(1.0+A*cos(k*x));
}

//Interesante pero no tan útil de implementar en 1D.
double einasto(double x, double v, double sx, double sv, double amplitude)
{
        //double x0 = 1.0; //EN vía lactea 20kpc.
        double gamma = 0.17;
        double in = pow(x/sx,gamma)-1.0;
        double exx = -2.0*in/gamma;
        //double eina = amplitude*exp(exx);
        //double exv = -v*v/(2*sv*sv);
        //double max = sqrt(2/PI)*exp(exv); //Para la velocidad se una una gaussiana
        return exx;
    
}

//Calcula la densidad, la velocidad macroscópica (u) y la energía libre (e) y los carga en los correspondientes arreglos.
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
                energy[i] += phase[i][j]*dv*pow(giveVel(j) - velocity[i], 2)/2.0;
			}
        velocity[i] = velocity[i] / density[i];
        energy[i] = energy[i] / density[i];
		mass += density[i]*dx;
		}
	return mass;
}

//Calcula el potencial (V) con el método de Fourier. Actualiza el arreglo pot.
void potencial()
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

        in[i] = giveDensity(i);
        inR[i] = -1.0;
        out[i] = 0;
        //in[i] = sin(2.0*PI*i*deltax);//Actualmente funciona para sin(x) confirmado. (se esperaba que la parte real de out fuera 0 y la imaginaria tuviera los picos, sin embargo solo el módulo cumple esto)
    }
    fftw_execute(pIda);

    //Guarda out en mem.
    for(i=0;i<Nx;i+=1){
        mem[i] = out[i];
    }

    pIda = fftw_plan_dft_1d(Nx, out, inR, FFTW_BACKWARD, FFTW_MEASURE); //Se debe usar el mismo plan sí o sí al parecer.
    
    //Devuelve carga a out Î(Chi).
    out[0] = -4*PI*G*mem[0];
    for(i=1;i<Nx;i+=1){
      out[i] = -4*PI*G*mem[i]/calcK2((double)i);
    //out[i] = mem[i]; //Descomentar esta línea para obtener la distribucion original.
    }
    fftw_execute(pIda);

    
    for(i=0;i<Nx;i+=1){
        pot[i] = creal(inR[i]/Nx);
    //    fprintf(oR, "%f\n",creal(inR[i])/Nx);
    //    fprintf(oI, "%f\n",cimag(inR[i])/Nx);
    }
}

//Calcula el k**2 de mis notas.
double calcK2(double j2)
{
    if((j2 == Nx/2.0)){
        return 1.0;
    }
    double k2= 2.0*sin(dx*PI*j2)/dx;
    return pow(k2,2);
}

//Método para evitar efectos misticos en la memoria.
double giveDensity(int l)
{
    double rta = density[l];
    return rta;
}

//Convierte unidades de la simulación a masas solares, metros, o segundos.
double convertir(double valor, int unidad )
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

//Deriva el potencial y carga la aceleración gravitacional en el arreglo acce.
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
            fprintf(output, "%f\n",convertir(acce[i], aKpc)/pow(convertir(1.0, aByear)*1000,2)); //Imprime en kpc / s^2;
			}
	fclose(output);
}

//Imprime el arreglo Pot
void printPot(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",pow(convertir(pot[i], aMetros)/convertir(1.0, aSegundos),2)/pot[i]);
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
//        if(i2 >= Nx){
//            printf("i = %d\n", i2);
//        }
        v = Vmin*1.0+dv*j2;
        
        //parte colisional.
//        v += collision(iin, jin, TAU);
        //
        x = v*dt;
        double di = x/dx*scale;
        di = (int) di;

        i2 = iin + di;
        i2 = mod(i2,Nx);
//	printf("%d\n",j2);
    return 0;
}

//Calcula un paso. Guarda una copia del phase actual en phaseOld. Actualiza phase. k,i son x. j,l son v.
void step()
{
	for(k = 0; k<Nx; k++){
		for(l= 0; l<Nv; l++){
			if(newij(k,l) ==0){
				//phaseOld[k][l] = phase[k][l]; TODO:No se usa
				phaseTemp[i2][j2] += phase[k][l];
                    //phaseTemp[i2][l] += collision(k,mod(l,1024),TAU);
                    //phaseTemp[k][l] -= collision(k,l,TAU);
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
				//phaseOld[k][l] = phase[k][l];
				//phaseTemp[i2][j2] += phase[k][l];
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
        double dj = v/dv;
        dj = (int)dj;
        j2 = jin;

        if(j2 < 0 || j2 >= Nv) return -1;
//        if(i2 >= Nx){
//            printf("i = %d\n", i2);
//        }
        v = giveVel(j2);
        
        //parte colisional.
//        v += collision(iin, jin, TAU);
        //
        x = v*dt*scale;
        double di = x/dx;
        di = (int) di;

        i2 = iin + di;
        i2 = mod(i2,Nx);
//	printf("%d\n",j2);
    return 0;
}

//Calcula la contribución colisional en phase[icol][jcol] con un Tau dado.
double collision(int icol, int jcol, double tau)
{
    if(TAU==0) return 0;
    double df = (feq(icol,jcol) - phaseOld[icol][jcol])/tau;    
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

















