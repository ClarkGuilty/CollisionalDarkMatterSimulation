/*
Javier Alejandro Acevedo Barroso

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
double *phaseTemp;
double *density;
double *pot;
double *accex;
double *accey;

//Variables recurrentes
int i;
int j;
int k;
int l;

//Para las nuevas posiciones.
int i2x;
int j2x;
int i2y;
int j2y;

//Para iterar.
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
double potencial(); 
double calcK2(double i2, double j2);
double convertir(double valor, int unidad);
void calAcce(); 
void printAcce(char *name); //TODO
double newij(int iinx, int jinx, int iiny int jiny); //TODO
void step(); //TODO
int mod(int p, int q);
void printPot(char *name);
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
    accex = malloc((sizeof(double)*Nx));
    accey = malloc((sizeof(double)*Nx));
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
        calAcce();
        
        

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
            //printf("%f\n", mem[in(k1,k2)]);
            //printf("%f %f %d %d\n",calcK2((double)k1,(double)k2), mem[in(k1,k2)],k1,k2);            
            out[in(k1,k2)] = -PI*G*mem[in(k1,k2)]*calcK2((double)k1,(double)k2);//Porque dx = dy y Nx = Ny.
        //out[in(k1,k2)] = mem[in(k1,k2)]; //Descomentar esta línea para obtener la distribucion original.
    
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

//Retorna 1/(sin2 + sin2).
double calcK2(double i2, double j2)
{
    if( ( (j2 == 0) || (j2 == Nx/2) )  && ( (i2 == 0) || (i2 == Nx/2) )  ){
        return 0;
    }
    double rta1= sin(dx*PI*j2);
    double rta2= sin(dx*PI*i2);
    return 1.0/(rta1*rta1+rta2*rta2);
}

//Retorna la densidad en (in1,in2).
double giveDensity(int in1, int in2)
{
    double rta = density[in(in1,in2)];
    return rta;
}

//Retorna el potencial en (in1,in2).
double givePot(int in1, int in2)
{
    double rta = pot[in(in1,in2)];
    return rta;
}

//Retorna aceleración en x en (in1,in2).
double giveAccex(int in1, int in2)
{
    double rta = accex[in(in1,in2)];
    return rta;
}

//Retorna aceleración en y en (in1,in2).
double giveAccey(int in1, int in2)
{
    double rta = accey[in(in1,in2)];
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
}

//Imprime el arreglo Acce.
void printAcce(char *name)
{
	FILE *outputx = fopen('x'+name, "w+");
    FILE *outputy = fopen('y'+name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(outputx, "%f\n",accex[i]);
            fprintf(outputy, "%f\n",accey[i]);
			}
	fclose(outputx);
    fclose(outputy);
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
/        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);
}


//Calcula la nueva posición de (x,y,vx,vy)=(iinx,iiny,jinx,jiny).
double newij(int iinx, int iiny,int jinx, int jiny)
{
        scale = 1.0; //1.0 es el valor estandar, modificarlo para mejorar visualización.
        //double x = Xmin*1.0+dx*iinx; //Inicialización
        double vx = giveAccex(iinx,iiny)*dt; //Cambio de velocidad.
        double djx = vx/dvx;                
        djx = (int)djx;                     //Cambio de vel en el tablero.
        j2x = jinx+djx;                     //Nuevo j.

        //double y = Ymin*1.0+dy*iiny; //Inicialización
        double vy = giveAccey(iiny,iiny)*dt;
        double djy = vy/dvy;
        djy = (int)djy;
        j2y = jiny+djy;
        
        if(j2x < 0 || j2x >= Nvx || j2y < 0 || j2y >= Nvy) return -1;
//        if(i2 >= Nx){
//            printf("i = %d\n", i2);
//        }
        vx = Vxmin*1.0+dvx*j2x;
        vy = Vymin*1.0+dvy*j2y;
        
        //parte colisional.
//        v += collision(iin, jin, TAU);
        //
        x = vx*dt;
        double dix = x/dx*scale;
        dix = (int) dix;
        
        y = vy*dt;
        double diy = y/dy*scale;
        diy = (int) diy;
        
        

        i2x = iinx + dix;
        i2x = mod(i2x,Nx);
        
        i2y = iiny + diy;
        i2y = mod(i2y,Ny);
//	printf("%d\n",j2);
    return 0;
}


//Calcula un paso. Guarda una copia del phase actual en phaseOld. Actualiza phase. 
void step()
{
    
    
    
    for(k1=0;k1<Nx;k1+=1) {
        for(k2=0;k2<Ny;k2+=1) {
            for(k3=0;k3<Nvx;k3+=1) {
                for(k4=0;k4<Nvy;k4+=1) {
                    if(newij(k1,k2,k3,k4) == 0){
                        phaseTemp[ind(i2x,i2y,j2x,j2y)] += phase[ind(k1,k2,k3,k4)];
                    }
                }
            }
        }
    }
    
    
    
    
    for(k1=0;k1<Nx;k1+=1) {
        for(k2=0;k2<Ny;k2+=1) {
            for(k3=0;k3<Nvx;k3+=1) {
                for(k4=0;k4<Nvy;k4+=1) {
                    phase[ind(k1,k2,k3,k4)] = phaseTemp[ind(k1,k2,k3,k4)];
                    phaseTemp[ind(k1,k2,k3,k4)] = 0;
                }
            }
        }
    }
    
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

















