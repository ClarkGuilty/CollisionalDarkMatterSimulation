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
//#define Xmin -1.0
//#define Xmax 1.0
//#define Ymin -1.0
//#define Ymax 1.0
#define Xmin 0
#define Xmax 2.0
#define Ymin 0
#define Ymax 2.0

#define Vxmin -1.0
#define Vxmax 1.0
#define Vymin -1.0
#define Vymax 1.0


//Tamaño del espacio.
#define tamano 128
#define Nx tamano
#define Ny tamano
#define Nvx tamano
#define Nvy tamano


//Constantes de unidades.
#define aMetros 18
#define aSegundos 14
#define aByear 4
#define aMasasSol 5




#define mParsecs 200e-3  //Cuántos megaparsecs equivalen a una unidad espacial.
#define solarMases 1e12 //Cuántas masas solares equivalen a una unidad de masa.
#define fracT0 3e-3     //Qué fracción de la edad del universo equivale a una unidad de tiempo
//#define G 0.00096 //G en estas unidades. Se calcula con calculations.py
#define G 1.0

#define radius 0.5
#define MASSo 10.0
#define scale 1.0 //1.0 es el valor estandar, modificarlo para mejorar visualización.



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
double *densitytheo;
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

double totalMass = 0;

fftw_complex *inE, *out, *inR, *mem, *out2;
fftw_plan pIda;
fftw_plan pVuelta;

double ix;
double iy;
double ivx;
double ivy;

double deltax;
double deltay;
double deltavx;
double deltavy;

double dix;
double diy;
double djx;
double djy;

double Lx = Xmax - Xmin;
double Ly = Ymax - Ymin;
double Lvx = Vxmax - Vxmin;
double Lvy = Vymax - Vymin;
double dx = (Xmax - Xmin)*1.0/(Nx);
double dy = (Ymax - Ymin)*1.0/(Ny);
double dvx = (Vxmax - Vxmin)*1.0/(Nvx);
double dvy = (Vymax - Vymin)*1.0/(Nvy);

double dt = 0.5;
int Nt = 15;

double totalPerdido;

double sigma = MASSo/(PI*radius*radius);
double lambda = MASSo/2;


FILE *constantes;
void printPhase(char *name); //TODO: debe replantearse lo que se va a imprimir.
void printPhaseX(char *name, int corteY, int corteVy);
void printPhaseY(char *name, int corteX, int corteVx);
double gaussD(double x,double y, double vx, double vy, double sr, double sv, double amplitude); 
double calDensity();
void printDensity(char *name);
void printConstant(char *name, double value);
double giveDensity(int in1,int in2);
double giveAccex(int in1, int in2);
double giveAccey(int in1, int in2);
double givePot(int in1, int in2);
double potencial(); 
double calcK2(double i2, double j2);
double convertir(double valor, int unidad);
void calAcce(); 
void printAccex(char *namex);
void printAccey(char *namex);
double newij(int iinx, int jinx, int iiny, int jiny); 
void step(); 
int mod(int p, int q);
void printPot(char *name);
int ind(int in1, int in2, int in3, int in4);
int in(int in1, int in2);
double gtheory(int in1, int in2,int eje);
double darX(int input);
double darY(int input);
double darVx(int input);
double darVy(int input);
double gmod(int in1, int in2);
void printGtheo(char *namex, char *namey);
double uniDisk(double x, double y);
double densidadTeorica(int inx, int iny);
void loadPot(double alpha, int n);
void printDensityTheo(char *name);
void loadDensityTheo();
double newPot(double alpha, int n, int iin, int jin);
double laplacePot(int iin, int jin,double alpha, int n);
double potencialTeorico(int inx, int iny);
double potencialTeorico2(int inx, int iny, int nx, int ny);
double potencial2(); 
double densidadTeorica2(int inx, int iny, int nx, int ny);
double calcK3(double i2, double j2);
double calcK4(double i2, double j2);

int main()
{
    dt = dt*dx*dy/dvx/dvy;
    phase = malloc((sizeof(double)*Nx*Ny*Nvx*Nvy));
    phaseTemp = malloc((sizeof(double)*Nx*Ny*Nvx*Nvy));
    if(phase == NULL){
        printf("phase es Null\n");   
        }
    density = malloc((sizeof(double)*Nx*Ny));
    densitytheo = malloc((sizeof(double)*Nx*Ny));
    accex = malloc((sizeof(double)*Nx*Ny));
    accey = malloc((sizeof(double)*Nx*Ny));
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
	double sr = 0.13;
    double sv = 0.13;
	double ampl = 0.0001;
    int nx = 10;
    int ny = 4;
        //printf("size of double %lu\n", sizeof(double));
        printf("%d %d %d %d\n", Nx,Ny,Nvx,Nvy);
        //phase[0][0][0][1] = 1;
        //double mass0 = calDensity();
        
        //printf("Masa = %f\n", mass0);
        
//        radius = 0.1;
        totalMass = 0;
        for(k1=0;k1<Nx;k1+=1) {
            x = Xmin*1.0+dx*k1;
            for(k2=0;k2<Ny;k2+=1) {
                y = Ymin*1.0+ dy*k2;
                
                //density[in(k1,k2)] = potencialTeorico2(k1,k2,3,20);
                density[in(k1,k2)] = densidadTeorica2(k1,k2,nx,ny);

                pot[in(k1,k2)] = potencialTeorico2(k1,k2,nx,ny);
                totalMass += density[in(k1,k2)]*dx*dy;// + 2.0/PI/Nx/Ny;
            }
            //printf("%d\n",k1);
        }
        
        
        
        //loadPot(4*PI*G*MASSo,2);
        
        printDensity("./datFiles/density0.dat");
        printPot("./datFiles/potential0.dat");
        potencial2();
        printAccex("./datFiles/fpot0.dat");
        printAccey("./datFiles/fpot1.dat");
        
        printPot("./datFiles/potential1.dat");
        

    fclose(constantes);
	//fclose(simInfo);
	return 0;

}


//Calcula el potencial (V) con el método de Fourier. Actualiza el arreglo pot.
double potencial()
{
    inE=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    inR=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    mem=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    pIda = fftw_plan_dft_2d(Nx, Ny, inE, out,FFTW_FORWARD, FFTW_MEASURE);

    //Cargar densidad en in:
    for(k1=0;k1<Nx;k1+=1){
        for(k2=0;k2<Ny;k2+=1){
        inE[in(k1,k2)] = giveDensity(k1,k2)- totalMass/((Xmax-Xmin)*(Ymax-Ymin));
        inR[in(k1,k2)] = -1.0;
        out[in(k1,k2)] = 0;
        }
    }

    
    fftw_execute(pIda);
    
    //Guarda out en mem. 
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2<Ny; k2+=1){
         mem[in(k1,k2)] = out[in(k1,k2)];   
        }
    }

    pIda = fftw_plan_dft_2d(Nx, Ny, out, inR,FFTW_BACKWARD, FFTW_MEASURE);
    //Se debe usar el mismo plan sí o sí al parecer.

    //Devuelve carga a out Î(Chi).
    
    //out[0] = -4*PI*G*mem[0] ;
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny   ;k2 += 1){
            //printf("%f\n", mem[in(k1,k2)]);
            //printf("%f %f %d %d\n",calcK2((double)k1,(double)k2), mem[in(k1,k2)],k1,k2);            
            out[in(k1,k2)] = -4*PI*G*(mem[in(k1,k2)])*calcK2((double)k1,(double)(k2));//Porque dx = dy y Nx = Ny.
            //printf("%f\n", mem[in(k1,k2)]);
            //out[in(k1,k2)] = -PI*G*mem[in(k1,k2)]/(creal(cpow(I*PI*dx*k1,2))+ creal(cpow(I*PI*dy*k2,2)));//Random plan b que parece fallar.
            //printf("%f\n", creal(cpow(I*PI*dx*k2,2)));
            
        //out[in(k1,k2)] = mem[in(k1,k2)]; //Descomentar esta línea para obtener la distribucion original.
        }
    }
    
    
        k1 = 0;
        printf("a0 es %f\n", creal(mem[in(0,0)]));
     //   out[in(0,0)] = 0;
        
        
    for(k2 = 1; k2<Ny;k2+=1){
      //  out[in(k1,k2)] = -PI*G*(mem[in(k1,k2)])*calcK2((double)k1,(double)(k2));//Porque dx = dy y Nx = Ny.
    }

//out[0] = PI*G*totalMass; //Posible alternativa para a0. Doesn't seems like it.
    fftw_execute(pIda);
    printf("%f\n",totalMass);


    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1)
        pot[in(k1,k2)] = creal(inR[in(k1,k2)]/Nx/Ny);

    }

 
 fftw_free(inE);
 fftw_free(out);
 fftw_free(inR);
 fftw_free(mem);
 fftw_free(out2);
 
}


//Calcula el potencial (V) con el método de Fourier. Actualiza el arreglo pot.
double potencial2()
{
    inE=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    inR=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    mem=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    pIda = fftw_plan_dft_2d(Nx, Ny, inE, out, FFTW_FORWARD, FFTW_MEASURE);

    //Cargar densidad en in:
    for(k1=0;k1<Nx;k1+=1){
        for(k2=0;k2<Ny;k2+=1){
        inE[in(k1,k2)] = giveDensity(k1,k2);//- totalMass/((Xmax-Xmin)*(Ymax-Ymin));
        inR[in(k1,k2)] = 0;
        out[in(k1,k2)] = 0;
        }
    }

    
    //fftw_execute(pIda);
    fftw_execute(pIda);
    
    //Guarda out en mem. 
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2<Ny; k2+=1){
         mem[in(k1,k2)] = out[in(k1,k2)];   
        }
    }

    pIda = fftw_plan_dft_2d(Nx, Ny, out, inR,FFTW_BACKWARD, FFTW_MEASURE);
    //Se debe usar el mismo plan sí o sí al parecer.

    //Devuelve carga a out Î(Chi).
    
        for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1){
        accex[in(k1,k2)] = creal(mem[in(k1,k2)]/Nx/Ny); //Voy a guardar en accex mi espacio de fourier
        accey[in(k1,k2)] = cimag(mem[in(k1,k2)]/Nx/Ny); //Voy a guardar en accex mi espacio de fourier
        }
    }

        printAccex("./datFiles/fdens0.dat");
        printAccey("./datFiles/fdens1.dat");

    
    //out[0] = -4*PI*G*mem[0] ;
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny   ;k2 += 1){
            //printf("%f\n", mem[in(k1,k2)]);
            //printf("%f %f %d %d\n",calcK2((double)k1,(double)k2), mem[in(k1,k2)],k1,k2);            
             out[in(k1,k2)] = -4.0*PI*1.0*(mem[in(k1,k2)])*calcK4((double)k1,(double)k2);//Porque dx = dy y Nx = Ny.
            //printf("%f\n", mem[in(k1,k2)]);
            //out[in(k1,k2)] = -PI*G*mem[in(k1,k2)]/(creal(cpow(I*PI*dx*k1,2))+ creal(cpow(I*PI*dy*k2,2)));//Random plan b que parece fallar.
            //printf("%f\n", creal(cpow(I*PI*dx*k2,2)));
            
        //out[in(k1,k2)] = mem[in(k1,k2)]; //Descomentar esta línea para obtener la distribucion original.
        }
    }
    

     for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1){
        accex[in(k1,k2)] = creal(out[in(k1,k2)]); //Voy a guardar en accex mi espacio de fourier
        accey[in(k1,k2)] = cimag(out[in(k1,k2)]); //Voy a guardar en accex mi espacio de fourier
        }
    }
    
        k1 = 0;
        printf("a0 es %f\n", creal(mem[in(0,0)]));
     //   out[in(0,0)] = 0;




//out[0] = PI*G*totalMass; //Posible alternativa para a0. Doesn't seems like it.
    fftw_execute(pIda);
    printf("Masa total %f\n",totalMass);


    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1)
        pot[in(k1,k2)] = creal(inR[in(k1,k2)]/Nx/Ny);

    }

 
 fftw_free(inE);
 fftw_free(out);
 fftw_free(inR);
 fftw_free(mem);
 fftw_free(out2);
 
}

//Retorna 1/(sin2 + sin2).
double calcK2(double i2, double j2)
{
    //if( ( (j2 == 0) || (j2 == Nx/2) )  && ( (i2 == 0) || (i2 == Nx/2) )  ){
    if( ( (j2 == 0) )  && ( (i2 == 0)   )  ){
        return 1;
    }
     if(i2<Nx/2-1){
         i2 = PI*i2;
     }
     if(j2<Nx/2-1){
         j2 = PI*j2;
     }
     if(i2>=Nx/2-1){
         i2 = PI*(i2-Nx);
     }
     if(j2>=Nx/2-1){
         j2 = PI*(j2-Ny);
     }
    //double rta1= sin(dx*PI*j2);
    //double rta2= sin(dx*PI*i2);
//     double rta1= -4*Nx*sin(PI*i2/Nx)/PI;
//     double rta2= -4*Nx*sin(PI*j2/Nx)/PI;

    //return 1.0/(pow(rta1,2)+pow(rta2,2));
    return -1.0/(pow(i2,2)+pow(j2,2));
}

double calcK3(double i2, double j2)
{
    //if( ( (j2 == 0) || (j2 == Nx/2) )  && ( (i2 == 0) || (i2 == Nx/2) )  ){
    if( ( (j2 == 0) || (j2==  Nx/2) )  && ( (i2 == 0) || (i2 == Nx/2)   )  ){
        return 0;
    }
     if(i2<Nx/2+1){
         i2 = i2;
     }
     if(j2<Nx/2+1){
         j2 = j2;
     }
     if(i2>=Nx/2+1){
         i2 = -i2;
     }
     if(j2>=Nx/2+1){
         j2 = -j2;
     }
    double rta1= 2*sin(dx*j2*PI)/dx;
    double rta2= 2*sin(dx*i2*PI)/dx;
//     double rta1= -4*Nx*sin(PI*i2/Nx)/PI;
//     double rta2= -4*Nx*sin(PI*j2/Nx)/PI;

    
    //printf(" (%.0f, %.0f) = %f\n",i2,j2, -1.0/(pow(rta1,2)+pow(rta2,2)));
    return 1.0/(pow(rta1,2)+pow(rta2,2));
    //return 1.0/(i2*i2+j2*j2);

}


double calcK4(double i2, double j2)
{
    if( ( (j2 == 0) || (j2==  Nx/2) )  && ( (i2 == 0) || (i2 == Nx/2)   )  ){
        return 0;
    }
     if(i2<Nx/2+1){
         i2 = PI*i2;
     }
     if(j2<Nx/2+1){
         j2 = PI*j2;
     }
     if(i2>=Nx/2+1){
         i2 = -PI*(Nx-i2);
     }
     if(j2>=Nx/2+1){
         j2 = -PI*(Nx-j2);
     }    
    return 1.0/(i2*i2+j2*j2);

}


double densidadTeorica2(int inx, int iny, int nx, int ny)
{
    
 return -(pow(PI*0.5*2.0*nx,2)*sin(PI*darX(inx)*0.5*nx*2.0)+pow(PI*0.5*ny*2.0,2)*sin(PI*darY(iny)*0.5*ny*2.0))/(4.0*PI);

}

//POtencial con sin(nx*pi*x/Nx) + sin(ny*pi*y/Ny)
double potencialTeorico2(int inx, int iny, int nx, int ny)
{
    
 return sin(PI*darX(inx)*0.5*nx*2)+sin(PI*darY(iny)*0.5*ny*2.0);
 //np.sin(np.pi*xin*nx*2.0/L)+np.sin(np.pi*yin*ny*2.0/L)

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

//Imprime el corte y = corteY , Vy = corteVy del espacio de fase (un plano 2d).
void printPhaseX(char *name, int corteY, int corteVy)
{
 
   	FILE *output = fopen(name, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Nvx+1;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", convertir(phase[ind(i,Ny/2,Nvx-j,Nvy/2)], aMasasSol)/convertir(1.0,aKpc)/(convertir(1.0,aKpc)*3.0857e+19)* convertir(1.0,aSegundos)); //Imprime en Masas solares /kpc / (km/s)
            fprintf(output,"%f ",phase[ind(i,corteY,Nvx-j,corteVy)]);
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);

    
}

//Imprime el corte x = corteX , Vx = corteVx del espacio de fase (un plano 2d).
void printPhaseY(char *name, int corteX, int corteVx)
{
 
   	FILE *output = fopen(name, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Nvx+1;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", convertir(phase[ind(i,Ny/2,Nvx-j,Nvy/2)], aMasasSol)/convertir(1.0,aKpc)/(convertir(1.0,aKpc)*3.0857e+19)* convertir(1.0,aSegundos)); //Imprime en Masas solares /kpc / (km/s)
            fprintf(output,"%f ",phase[ind(corteX,i,corteVx,Nvy-j)]);
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
		for(j=0;j<Ny;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", phase[i][Nv-j]);
                    //fprintf(output,"%f ", giveDensity(i,mod(Ny-j,Ny)));  //The correct way to print. Use .T in numpy.
            fprintf(output,"%f ", giveDensity(i,j));  //Al parecer esta SÍ es la real.
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);

}

//Imprime el arreglo density con el String name como nombre.
void printAccex(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", phase[i][Nv-j]);
                    //fprintf(output,"%f ", giveAccex(i,mod(Ny-j,Ny)));  //The correct way to print. Use .T in numpy.
            fprintf(output,"%f ", giveAccex(i,j));
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);

}

//Imprime el arreglo density con el String name como nombre.
void printAccey(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", phase[i][Nv-j]);
                    //fprintf(output,"%f ", giveAccey(i,mod(Ny-j,Ny)));  //The correct way to print. Use .T in numpy.
                    fprintf(output,"%f ", giveAccey(i,j));
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

double uniDisk(double x, double y)
{
 if(x*x+y*y < radius*radius)
 {
  return sigma;
     
 }
 
 return 0;
    
}



//Calcula la densidad. Actualiza el arreglo density
double calDensity()
{
    totalMass = 0;
        for(k1=0; k1<Nx;k1+=1){
            for(k2=0; k2< Ny; k2+=1){
                density[in(k1,k2)] = 0;
                for(k3=0; k3< Nvx; k3+=1){
                    for(k4=0; k4< Nvy; k4+=1){
                        //density[in(k1,k2)] += phase[k1][k2][k3][k4]*dvx*dvy;
                        density[in(k1,k2)] = giveDensity(k1,k2)+ phase[ind(k1,k2,k3,k4)];
                    }
                }
                totalMass += density[in(k1,k2)];
            }
        }
        return totalMass;
}




//double calcK2(double i2, double j2)
//{
//     if(i2<=Nx/2){
//         i2 = PI*i2;
//     }
//     if(j2<=Nx/2){
//         j2 = PI*j2;
//     }
//     if(i2>Nx/2){
//         i2 = PI*(i2-Nx);
//     }
//     if(j2>Nx/2){
//         j2 = PI*(j2-Ny);
//     }
// 
//     return 1.0/(pow(i2,2) + pow(j2,2));
//     
//}

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
        for(k2 = 0; k2<Ny ; k2 +=1){
    accex[in(k1,k2)] =  (-givePot(mod(k1-2,Nx),k2) + 8*givePot(mod(k1-1,Nx),k2)-8*givePot(mod(k1+1,Nx),k2)+givePot(mod(k1+2,Nx),k2))/(12*dx);
    accey[in(k1,k2)] =  (-givePot(k1,mod(k2-2,Ny)) + 8*givePot(k1,mod(k2-1,Ny))-8*givePot(k1,mod(k2+1,Ny))+givePot(k1,mod(k2+2,Ny)))/(12*dy);
            //printf("%d %d\n" ,k1, mod(k2-2,Ny));
    
    }
}
}


void printGtheo(char *namex, char *namey)
{
	FILE *outputx = fopen(namex, "w+");
    FILE *outputy = fopen(namey, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Ny+1;j+=1){ 
          
            fprintf(outputx, "%f ",gtheory(i,Ny-j,0));
            fprintf(outputy, "%f ",gtheory(i,Ny-j,1));
        
        }
		fprintf(outputx,"\n");
        fprintf(outputy,"\n");
		//printf("%d\n", i);
        }
	fclose(outputx);
    fclose(outputy);
}

//Imprime el arreglo Pot
void printPot(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
          //      printf("ignorarPrimero\n");
			//fprintf(output,"%f ", phase[i][Nv-j]);
                    fprintf(output,"%f ", givePot(i,mod(Ny-j,Ny)));
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);
}


//Calcula la nueva posición de (x,y,vx,vy)=(iinx,iiny,jinx,jiny).
double newij(int iinx, int iiny,int jinx, int jiny)
{
        
        //double x = Xmin*1.0+dx*iinx; //Inicialización
        deltavx = giveAccex(iinx,iiny)*dt; //Cambio de velocidad.
        djx = deltavx/dvx;                
        djx = (int)djx;                     //Cambio de vel en el tablero.
        j2x = jinx+djx;                     //Nuevo j.

        //double y = Ymin*1.0+dy*iiny; //Inicialización
        deltavy = giveAccey(iinx,iiny)*dt;
        djy = deltavy/dvy;
        djy = (int)djy;
        j2y = jiny+djy;
        
        if(j2x < 0 || j2x >= Nvx || j2y < 0 || j2y >= Nvy)
        {
            totalPerdido += phase[ind(iinx,iiny,jinx,jiny)];
            return -1;
            
        }
//        if(i2 >= Nx){
//            printf("i = %d\n", i2);
//        }
        deltavx = Vxmin*1.0+dvx*j2x;
        deltavy = Vymin*1.0+dvy*j2y;
        
        //parte colisional.
//        v += collision(iin, jin, TAU);
        //
        deltax = deltavx*dt;
        dix = deltax/dx*scale;
        dix = (int) dix;
        
        deltay = deltavy*dt;
        diy = deltay/dy*scale;
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



//Da el valor de la gravedad para una distribución de disco en (in1,in2) en la cordenada dada en eje (x = 0, y = 1).
double gtheory(int in1, int in2, int eje)
{
    double angle = 0;
    if(darX(in1) == 0)
    {
        if(eje == 1)
        {
         return -gmod(in1,in2)*fabs(darY(in2))/darY(in2);   
        }
        return 0;
        
    }
    if(eje == 1)
    {
     return -gmod(in1,in2)*sin(   atan(darY(in2)/darX(in1)) )* fabs(darX(in1))/darX(in1); //Modulo del vector * sin (theta), la dirección de y va incluida en sin.  atan va de -pi/2 a pi/2. 
       // return 0;
    }
    
    return -1*gmod(in1,in2)*cos( fabs( atan(darY(in2)/darX(in1)) ) )* fabs(darX(in1))/darX(in1);//Modulo del vector * cos(theta) * en dirección de x
    //return 2*PI*lambda*G;
    
    //return 0;
}


double gmod(int in1, int in2)
{
    
    if(darX(in1)*darX(in1)+darY(in2)*darY(in2)< radius*radius){
        
        return 2*G*sigma*PI*sqrt( darX(in1)*darX(in1) +darY(in2)*darY(in2));
    }
 return 2*G*MASSo/sqrt(darX(in1)*darX(in1)+darY(in2)*darY(in2));
}


double darX(int input)
{
 return Xmin*1.0+dx*input;   
}
double darY(int input)
{
 return Ymin*1.0+dy*input;   
}
double darVx(int input)
{
 return Vxmin*1.0+dvx*input;   
}
double darVy(int input)
{
 return Vymin*1.0+dvy*input;   
}



double newPot(double alpha, int n, int iin, int jin)
{
                    ix = darX(iin);
                    iy = darY(jin);
                    return -alpha * pow(sqrt(pow(ix,2)+pow(iy,2)),n);
}

//Calcula el laplaciano del potencial.
double laplacePot(int iin, int jin,double alpha, int n)
{
    //return (givePot(mod(iin+1,Nx),jin) + givePot(mod(iin-1,Nx),jin)+givePot(iin,mod(jin+1,Ny))+givePot(iin,mod(jin-1,Ny))-4*givePot(iin,jin))/(dy*dy); //dx = dy
//        printf("%d %d %f %f %f %f %f %f\n",iin,jin, givePot(mod(iin+1,Nx),mod(jin+1,Ny)),givePot(mod(iin-1,Nx),mod(jin+1,Ny)),givePot(mod(iin+1,Nx),mod(jin-1,Ny)),givePot(mod(iin-1,Nx),mod(jin-1,Ny)),givePot(iin,jin),    -(givePot(mod(iin+1,Nx),mod(jin+1,Ny)) +givePot(mod(iin-1,Nx),mod(jin+1,Ny))+givePot(mod(iin+1,Nx),mod(jin-1,Ny))+givePot(mod(iin-1,Nx),mod(jin-1,Ny))-4*givePot(iin,jin))    );
    //return -(givePot(mod(iin+1,Nx),mod(jin+1,Ny)) +givePot(mod(iin-1,Nx),mod(jin+1,Ny))+givePot(mod(iin+1,Nx),mod(jin-1,Ny))+givePot(mod(iin-1,Nx),mod(jin-1,Ny))-4*givePot(iin,jin))/(dy*dy); //dx = dy
    
return -(newPot(alpha, n,mod(iin,Nx),mod(jin-1,Ny)) +newPot(alpha, n,mod(iin,Nx),mod(jin+1,Ny))+newPot(alpha, n,mod(iin+1,Nx),mod(jin,Ny))+newPot(alpha, n,mod(iin-1,Nx),mod(jin,Ny))-4*newPot(alpha, n,iin,jin))/(dy*dy);
    
//return -(newPot(alpha, n,mod(iin+1,Nx),mod(jin+1,Ny)) +newPot(alpha, n,mod(iin-1,Nx),mod(jin+1,Ny))+newPot(alpha, n,mod(iin+1,Nx),mod(jin-1,Ny))+newPot(alpha, n,mod(iin-1,Nx),mod(jin-1,Ny))-4*newPot(alpha, n,iin,jin))/(dy*dy);
    

}

//Imprime el arreglo density con el String name como nombre.
void printDensityTheo(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
                  //fprintf(output,"%f ", densitytheo[in(i,Ny-j)]);
            //fprintf(output,"%f ", laplacePot(i,Ny-j));
            fprintf(output,"%f ", densitytheo[in(i,j)]);
        //printf("Error MesagenoIgno\n");
        }
		fprintf(output,"\n");
		//printf("%d\n", i);
			}

	fclose(output);

}
void loadDensityTheo(double alpha, int n)
{
	for(i=1;i<Nx;i+=1) {
		for(j=1;j<Ny;j+=1){ 
                  //fprintf(output,"%f ", densitytheo[in(i,Ny-j)]);
            //fprintf(output,"%f ", laplacePot(i,Ny-j));
            //fprintf(output,"%f ", laplacePot(i,j)/(4*PI*G));
            densitytheo[in(mod(i,Nx),mod(j,Ny))] = laplacePot(i,j,alpha, n);// /(4*PI*G);
            density[in(mod(i,Nx),mod(j,Ny))] = laplacePot(i,j,alpha, n);// /(4*PI*G);
        //printf("Error MesagenoIgno\n");
        }
		//printf("%d\n", i);
			}
			
    double valor = laplacePot(6,6,alpha, n);
    //Por las condiciones periódicas, esto se debe hacer manualmente.
    for(j = 0; j<Ny;j+=1){
        densitytheo[in(0,j)] = valor;// /(4*PI*G);
            density[in(0,j)] = valor;// /(4*PI*G);
    }
    
        for(i = 0; i<Ny;i+=1){
        densitytheo[in(i,0)] = valor;// /(4*PI*G);
            density[in(i,0)] = valor;// /(4*PI*G);
    }
    totalMass = laplacePot(6,6,alpha, n)*Nx*Ny;

}





