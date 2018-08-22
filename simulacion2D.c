/*
Javier Alejandro Acevedo Barroso
Primer Bosquejo. 2D con método de fourier.
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
#define Ymin -1.0
#define Ymax 1.0
#define Vxmin -1.0
#define Vxmax 1.0
#define Vymin -1.0
#define Vymax 1.0


//Tamaño del espacio.
#define Nx 128
#define Ny 128
#define Nvx 128
#define Nvy 128
#define Nv 1000

//Constantes de unidades.
#define aMetros 18
#define aSegundos 14
#define aByear 4
#define aMasasSol 5




#define mParsecs 20e-3  //Cuántos megaparsecs equivalen a una unidad espacial.
#define solarMases 1e11 //Cuántas masas solares equivalen a una unidad de masa.
#define fracT0 1e-3     //Qué fracción de la edad del universo equivale a una unidad de tiempo
#define G 0.010661906775769971 //G en estas unidades. Se calcula con sPlots.py

//Unidades funcionales para clusters galácticos.
//#define mParsecs 5
//#define solarMases 1e15
//#define fracT0 2e-1
//#define G 0.2729448134597113


//Arreglos
//static double phase[Nx][Ny][Nvx][Nvy] = {0};
double *phase;
double *density;
double *pot;
double *accex;
double *accey;

//Variables recurrentes
int i;
int j;
int k;
int l;
int i2;
int j2;
int k1;
int k2;
int k3;
int k4;

double Lx = Xmax - Xmin;
double Ly = Ymax - Ymin;
double Lvx = Vxmax - Vxmin;
double Lvy = Vymax - Vymin;
double dx = (Xmax - Xmin)*1.0/Nx;
double dy = (Ymax - Ymin)*1.0/Nx;
double dvx = (Vxmax - Vxmin)*1.0/Nvx;
double dvy = (Vymax - Vymin)*1.0/Nvy;

double dt = 0.5;
int Nt = 5;
FILE *constantes;
void printPhase(char *name); //TODO: debe replantearse lo que se va a imprimir.
double gaussD(double x,double y, double vx, double vy, double sr, double sv, double amplitude); 
double calDensity();
void printDensity(char *name);
void printConstant(char *name, double value);
double giveDensity(int in1,int in2);
double giveAccex(int in1, int in2);
double giveAccey(int in1, int in2);
double potencial(); //super TODO.
double calcK2(double j2); //TODO: estudiar la aproximación ahora bidimensionalmente.
double convertir(double valor, int unidad);
void calAcce(); //TODO: implementar derivada bidimensional.
void printAcce(char *name); //TODO
double newij(int iin, int jin); //TODO
void step(); //TODO
int mod(int p, int q);
void printPot(char *name); //TODO
double einasto(double x, double v, double sx, double sv, double amplitude); 
int ind(int in1, int in2, int in3, int in4);
int in(int in1, int in2);


int main()
{
    dt = dt*dx*dy/dvx/dvy;
    phase = malloc((sizeof(double)*Nx*Ny*Nvx*Nvy));
    if(phase == NULL){
        printf("phase es Null\n");   
        }
    density = malloc((sizeof(double)*Nx*Ny));
    acce = malloc((sizeof(double)*Nx));
    pot = malloc((sizeof(double)*Nx*Ny));

	constantes = fopen("constants.dat","w+");
	printConstant("Xmin",Xmin);
	printConstant("Ymin",Ymin);
	printConstant("Xmax",Xmax);
        printConstant("Ymax",Ymax);
	printConstant("Vxmin",Vxmin);
        printConstant("Vymin",Vymin);
	printConstant("Vxmax",Vxmax);
        printConstant("Vymax",Vymax);
	printConstant("Nx",Nx);
        printConstant("Ny",Ny);
	printConstant("Nvx",Nvx);
	printConstant("Nvy",Nvy);
	printConstant("Nt", Nt);
	double x;
	double vx;
        double y;
	double vy;
	double sr = 0.1;
        double sv = 0.1;
	double ampl = 10;
        //printf("size of double %lu\n", sizeof(double));
        printf("%d %d %d %d\n", Nx,Ny,Nvx,Nvy);
        //phase[0][0][0][1] = 1;
        //TODO Inicializar un espacio de fase 4D.
	for(k1=0;k1<Nx;k1+=1) {
            x = Xmin*1.0+dx*k1;
            for(k2=0;k2<Ny;k2+=1) {
                y = Ymin*1.0+ dy*k2;
                density[in(k1,k2)] = 0;
                for(k3=0;k3<Nvx;k3+=1) {
                    vx = Vxmin*1.0+ dvx*k3;
                    for(k4=0;k4<Nvy;k4+=1) {
                        vy = Vymin*1-0 + dvy*k4;
                        //printf("indices: %d %d %d %d\n", k1,k2,k3,k4);

                        //phase[k1][k2][k3][k4] = gaussD(x,y,vx,vy,sr,sv,ampl);
                        phase[ind(k1,k2,k3,k4)] = gaussD(x,y,vx,vy,sr,sv,ampl);
                        //printf("(%d,%d,%d,%d) %f\n", k1,k2,k3,k4,phase[ind(k1,k2,k3,k4)]);
                        //phaseTemp[ind(k1,k2,k3,k4)] = 0;
                        
                    }
                }
            }
            //printf("%d\n",k1);
        }
        printf("Masa = %f\n", calDensity());
        //printf("en 64 = %f\n", phase[ind(64,64,64,64)] );
        printDensity("Density0.dat");
        potencial();
        printPot("potAft.dat");
        

			fclose(constantes);
	return 0;

}



//Imprime el espacio de fase con el String name como nombre.
void printPhase(char *name)//no
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Ny+1;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", phase[i][Nv-j]);
                    fprintf(output,"%f ", giveDensity(i,Ny-j));
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
		for(j=1;j<Ny+1;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", phase[i][Nv-j]);
                    fprintf(output,"%f ", giveDensity(i,Ny-j));
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
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
double gaussD(double x,double y, double vx, double vy, double sr, double sv, double amplitude)
{
	//double ex = -x*x/(sr)-y*y/(sr)-vx*vx/(sv)-vy*vy/(sv);
    	double ex = -x*x/(2.0*sr*sr)-y*y/(2.0*sr*sr)-vx*vx/(2.0*sv*sv)-vy*vy/(2.0*sv*sv);
        
//	double ex = -x*x/(sx*sx)-v*v/(sv*sv);
        //printf("%f\n",ex);
	//return amplitude*exp(ex)/(2*PI*sx*sv);
	//return amplitude*exp(-sqrt(fabs(ex)));
        return amplitude*exp(ex);

}


//Interesante pero no tan útil de implementar en 1D.
double einasto(double x, double v, double sx, double sv, double amplitude)
{
        double x0 = 1.0; //EN vía lactea 20kpc.
        double gamma = 0.17;
        double in = pow(x/sx,gamma)-1.0;
        double exx = -2.0*in/gamma;
        double eina = amplitude*exp(exx);
        double exv = -v*v/(2*sv*sv);
        double max = sqrt(2/PI)*exp(exv); //Para la velocidad se una una gaussiana
        return exx;
    
}

//Calcula la densidad. Actualiza el arreglo density
double calDensity()
{
	double mass = 0;
        for(k1=0; k1<Nx;k1+=1){
            for(k2=0; k2< Ny; k2+=1){
                density[in(k1,k2)] = 0;
                for(k3=0; k3< Nvx; k3+=1){
                    for(k4=0; k4< Nvy; k4+=1){
                        //density[in(k1,k2)] += phase[k1][k2][k3][k4]*dvx*dvy;
                        density[in(k1,k2)] = giveDensity(k1,k2)+ phase[ind(k1,k2,k3,k4)];
                    }
                }
                mass += density[in(k1,k2)];
            }
        }
        return mass;
}

//Calcula el potencial (V) con el método de Fourier. Actualiza el arreglo pot.
double potencial()
{
//    double * densityIN= malloc(sizeof(double)*Nx);
//    for(i = 0;i<Nx;i+=1){
//        densityIN[i] = giveDensity(i);
//    }


    fftw_complex *inE, *out, *inR, *mem, *out2;
    inE=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    inR=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    mem=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);


    //fftw_complex *out2;
   // out2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    //double *in2;
    //in2 = (double*) malloc((sizeof(double)*Nx));

    fftw_plan pIda, pVuelta;
    pIda = fftw_plan_dft_2d(Nx, Ny, inE, out,FFTW_FORWARD, FFTW_MEASURE);

    //fftw_plan pIda2, pVuelta2;
    //pIda2 = fftw_plan_dft_r2c_1d(Nx, in2, out2, FFTW_MEASURE);

    //FILE *input = fopen("inF.dat", "w+");
    //FILE *output0 = fopen("outF0.dat", "w+");
    //FILE *output1 = fopen("outF1.dat", "w+");
    //FILE *oR = fopen("oR.dat", "w+");
    //FILE *oI = fopen("oI.dat", "w+");


    //Cargar densidad en in:
    for(k1=0;k1<Nx;k1+=1){
        for(k2=0;k2<Ny;k2+=1){
        inE[in(k1,k2)] = giveDensity(k1,k2);
       //in[i] = sin(2.0*PI*i*deltax);//Actualmente funciona para sin(x) confirmado. (se esperaba que la parte real de out fuera 0 y la imaginaria tuviera los picos, sin embargo solo el módulo cumple esto)
//        fprintf(input, "%f\n",creal(inE[i]));
        }
    }

    fftw_execute(pIda);
//    fftw_execute(pIda2);

    //Guarda out en mem. 
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2<Ny; k2+=1){
            
         mem[in(k1,k2)] = out[in(k1,k2)];   
            
            
            
            
        }
     //   fprintf(output0, "%f\n",creal(out[i]));
     //   fprintf(output1, "%f\n",cimag(out[i]));


//        mem[i] = out2[i];
//        fprintf(output0, "%f\n",creal(out2[i]));
//        fprintf(output1, "%f\n",cimag(out2[i]));

    }

    //fftw_execute(pIda);
    pIda = fftw_plan_dft_2d(Nx, Ny, out, inR,FFTW_BACKWARD, FFTW_MEASURE);
    //Se debe usar el mismo plan sí o sí al parecer.
    double memDo;

    //Devuelve carga a out Î(Chi).
    out[0] = mem[0];
    //out[0] = -4*PI*G*mem[0] ;
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1){
        out[in(k1,k2)] = -4*PI*G*mem[in(k1,k2)]/(calcK2((double)k1)+calcK2((double)k2));//Porque dx = dy y Nx = Ny.
        out[in(k1,k2)] = mem[in(k1,k2)]; //Descomentar esta línea para obtener la distribucion original.
    
        }

    }
    
    fftw_execute(pIda);


    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1)
        pot[in(k1,k2)] = creal(inR[in(k1,k2)]/Nx/Ny);
        //fprintf(oR, "%f\n",creal(inR[i])/(Nx*Ny));
        //fprintf(oI, "%f\n",cimag(inR[i])/Nx/Ny);
    }

 //   fclose(input);
 //   fclose(output0);
 //   fclose(output1);
 //   fclose(oR);
 //   fclose(oI);

    //fftw_destroy_plan(pVuelta);
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
double giveDensity(int in1, int in2)
{
    double rta = density[in(in1,in2)];
    return rta;
}

//Ahora para el potencial
double givePot(int in1, int in2)
{
    double rta = pot[in(in1,in2)];
    return rta;
}

//Retorna aceleración.
double giveAccex(int in1, int in2)
{
    double rta = density[accex(in1,in2)];
    return rta;
}

//
double giveAccey(int in1, int in2)
{
    double rta = density[accey(in1,in2)];
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
    if(unidad == aMetros){
        return valor * conx0* mParsecs;
    }
    if( unidad == aByear){
        return valor*13.772*fracT0;
    }
    if( unidad == aSegundos){
        return valor*cont0s*fracT0;
    }
}

//Deriva el potencial y carga la aceleración en el arreglo acce.
void calAcce()
{
    for(k1 = 0; k1<Nx ; k1 +=1){
        for(k2 = 0; k2<Nx ; k2 +=1){
    accex[in(k1,k2)] =  (-pot[in(mod(k1-2,Nx),k2)] + 8*pot[in(mod(k1-1,Nx),k2)]-8*pot[in(mod(k1+1,Nx),k2)]+pot[in(mod(k1+2,Nx),k2)])/(12*dx);
    accey[in(k1,k2)] =  (-pot[in(k1,mod(k2-2,Ny))] + 8*pot[in(k1,mod(k2-1,Ny))]-8*pot[in(k1,mod(k2+1,Ny))]+pot[in(k1,mod(k2+2,Ny))])/(12*dy);
    }
}

//Imprime el arreglo Acce.
void printAcce(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",acce[i]);
			}
	fclose(output);
}

//Imprime el arreglo Pot
void printPot(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Ny+1;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", phase[i][Nv-j]);
                    fprintf(output,"%f ", givePot(i,Ny-j));
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);
}


//Calcula la posición del elemento (in1,in2,in3,in4) del espacio de fase (x,y,vx,vy).
int ind(int in1, int in2, int in3, int in4)
{
    return in4 + Nvy*(in3 + Nvx*(in2 + Ny*in1));
    
}

//Calcula la posición del elemento (in1, in2) en los arreglos bidimensionales.
int in(int in1, int in2)
{
    return in2 + in1*Ny;
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

















