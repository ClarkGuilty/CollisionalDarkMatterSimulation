/*
Javier Alejandro Acevedo Barroso
Simulación del espacio de fase de un fluido tridimensional de Materia Oscura colisional con un tiempo de relajación TAU fijable.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>

//Constantes de la simulación.
#define PI 3.14159265359

//Valores límites para la posición y velocidad.
#define Xmin -1.0
#define Xmax 1.0
#define Ymin -1.0
#define Ymax 1.0
#define Zmin -1.0
#define Zmax 1.0

#define Vxmin -1.0
#define Vxmax 1.0
#define Vymin -1.0
#define Vymax 1.0
#define Vzmin -1.0
#define Vzmax 1.0


//Tamaño del espacio.
#define tamano 32
#define Nx tamano
#define Ny tamano
#define Nz tamano
#define Nvx tamano
#define Nvy tamano
#define Nvz tamano


//Constantes de unidades.
#define aMetros 18
#define aSegundos 14
#define aByear 4
#define aMasasSol 5




#define mParsecs 20e-3  //Cuántos megaparsecs equivalen a una unidad espacial.
#define solarMases 1e12 //Cuántas masas solares equivalen a una unidad de masa.
#define fracT0 3e-3     //Qué fracción de la edad del universo equivale a una unidad de tiempo
//#define G 0.00096 //G en estas unidades. Se calcula con calculations.py
#define G 0.959572

#define radius 0.5
#define MASSo 10.0
#define scale 1.0 //1.0 es el valor estandar, modificarlo para mejorar visualización.

#define TAU 0

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
double *velocityx;
double *velocityy;
double *velocityz;
double *energy;
double *densitytheo;
double *pot;
double *acce;
double *acceTheo;


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
int i2z;
int j2z;

//Para iterar.
int k1;
int k2;
int k3;
int k4;
int k5;
int k6;

//Para las funciones.
double totalMass = 0;
double totalPerdido;

fftw_complex *inE, *out, *inR, *mem, *out2;
fftw_plan pIda;
fftw_plan pVuelta;


double x;
double vx;
double y;
double vy;
double z;
double vz;

double ix;
double iy;
double iz;
double ivx;
double ivy;
double ivz;

double deltax;
double deltay;
double deltaz;
double deltavx;
double deltavy;
double deltavz;

double dix;
double diy;
double diz;
double djx;
double djy;
double djz;

double Lx = Xmax - Xmin;
double Ly = Ymax - Ymin;
double Lz = Zmax - Zmin;
double Lvx = Vxmax - Vxmin;
double Lvy = Vymax - Vymin;
double Lvz = Vzmax - Vzmin;
double dx = (Xmax - Xmin)*1.0/(Nx);
double dy = (Ymax - Ymin)*1.0/(Ny);
double dz = (Zmax - Zmin)*1.0/(Nz);
double dvx = (Vxmax - Vxmin)*1.0/(Nvx);
double dvy = (Vymax - Vymin)*1.0/(Nvy);
double dvz = (Vzmax - Vzmin)*1.0/(Nvz);

double dt = 0.5;
int Nt = 25;

double sigma = MASSo/(PI*radius*radius);
double lambda = MASSo/2;


FILE *constantes;
void printPhase(char *name); //TODO: debe replantearse lo que se va a imprimir. Se hará a medida que se vaya requiriendo.


void printConstant(char *name, double value);

void printDensityXY(char *name, int corteZ);
void printDensityYZ(char *name, int corteX);
void printDensityXZ(char *name, int corteY);

void printPotXY(char *name, int corteZ);
void printPotXZ(char *name, int corteY);
void printPotYZ(char *name, int corteX);

void printAcceXY(char *name, int corteZ, int xyz);
void printAcceYZ(char *name, int corteX, int xyz);
void printAcceXZ(char *name, int corteY, int xyz);
void printAcceXY2(char *name, int corteZ, int xyz);

void printThr(void (*f)(char*, int), char *name);
void printThr2(void (*f)(char*, int, int), char *name, int xyz);

double gaussD(double x,double y, double z, double vx, double vy, double vz, double sr, double sv, double amplitude); 
double gaussD2(double x,double y, double z, double sr, double sv, double amplitude);
double gaussD2Dens(double x,double y, double z, double sr, double sv, double amplitude);
double calDensity();
double giveDensity(int in1,int in2, int in3);
double giveAcce(int in1, int in2, int in3, int xyz);
double givePot(int in1, int in2, int in3);
double potencial();
double calcK2(double i2, double j2, double k2);
double convertir(double valor, int unidad);
void calAcce();
double newij(int iinx, int iiny, int iinz, int jinx, int jiny, int jinz);
void step();
int mod(int p, int q);
int ind(int in1, int in2, int in3, int in4, int in5, int in6);
int ina(int in1, int in2, int in3, int in4);
int in(int in1, int in2, int in3);
double darX(int input);
double darY(int input);
double darZ(int input);
double darVx(int input);
double darVy(int input);
double darVz(int input);
double densidadTeorica(int inx, int iny);
double potencialTeorico2(int inx, int iny, int inz, int nx, int ny, int nz);
double densidadTeorica2(int inx, int iny, int inz, int nx, int ny, int nz);
double calcK3(double i2, double j2, double k2);
double calcK4(double i2, double j2);
double darAcceTheo(double x, double y, double z, int xyz);
void collisionStep();
double newijCol(int iinx, int iiny, int iinz, int jinx, int jiny, int jinz);
double feq(int iposx, int iposy, int iposz, int jvelx, int jvely, int jvelz);
double collision(int icolx, int icoly, int icolz,  int jcolx, int jcoly, int jcolz, double tau);


double sr = 0.1;
double sv = 0.1;
double ampl = 1;


int main()
{
    dt = dt*dx*dy*dz/dvx/dvy/dvz;
    //phase = malloc((sizeof(double)*Nx*Ny*Nz*Nvx*Nvy*Nz));
    phase = malloc((sizeof(double)*Nx*Ny*Nz*Nvx*Nvy*Nvz));
    phaseTemp = malloc((sizeof(double)*Nx*Ny*Nz*Nvx*Nvy*Nvz));
    //phaseTemp = malloc((sizeof(double)*Nx*Ny*Nz*Nvx*Nvy*Nvz));
    if(phase == NULL){
        printf("phase es Null\n");   
        }
    density = malloc((sizeof(double)*Nx*Ny*Nz));
    densitytheo = malloc((sizeof(double)*Nx*Ny*Nz));
    acce = malloc((sizeof(double)*Nx*Ny*Nz*3));
    acceTheo = malloc((sizeof(double)*Nx*Ny*Nz*3));
    velocityx = malloc((sizeof(double)*Nx*Ny*Nz));
    velocityy = malloc((sizeof(double)*Nx*Ny*Nz));
    velocityz = malloc((sizeof(double)*Nx*Ny*Nz));
    energy = malloc((sizeof(double)*Nx*Ny*Nz));


    pot = malloc((sizeof(double)*Nx*Ny*Nz));

	constantes = fopen("./datFiles/constants.dat","w+");
	printConstant("Xmin",Xmin);
	printConstant("Ymin",Ymin);
    printConstant("Zmin",Zmin);
	printConstant("Xmax",Xmax);
    printConstant("Ymax",Ymax);
    printConstant("Zmax",Zmax);
	printConstant("Vxmin",Vxmin);
    printConstant("Vymin",Vymin);
    printConstant("Vzmin",Vzmin);
	printConstant("Vxmax",Vxmax);
    printConstant("Vymax",Vymax);
    printConstant("Vzmax",Vzmax);
	printConstant("Nx",Nx);
    printConstant("Ny",Ny);
    printConstant("Nz",Nz);
	printConstant("Nvx",Nvx);
	printConstant("Nvy",Nvy);
    printConstant("Nvz",Nvz);
	printConstant("Nt", Nt);
        printf("%d %d %d %d %d %d\n", Nx,Ny,Nz,Nvx,Nvy,Nvz);

    sr = 0.1;
    sv = 0.1;
    ampl = 1;
        
    for(k1=0;k1<Nx;k1+=1) {
        x =  Xmin*1.0+dx*k1;
            for(k2=0;k2<Ny;k2+=1) {
                y = Ymin*1.0+ dy*k2;
                for(k3=0;k3<Nz;k3+=1) {
                    z= darZ(k3);
                    density[in(k1,k2,k3)] = 0;
                    for(k4=0;k4<Nvx;k4+=1) {
                        vx = darVx(k4);
                        for(k5=0;k5<Nvy;k5+=1) {
                            vy = darVy(k5);
                            for(k6=0;k6<Nvz;k6+=1){
                                vz = darVz(k6);
                        //printf("indices: %d %d %d %d\n", k1,k2,k3,k4);

                        //phase[k1][k2][k3][k4] = gaussD(x,y,vx,vy,sr,sv,ampl);
                                phase[ind(k1,k2,k3,k4,k5,k6)] = gaussD(x,y,z,vx,vy,vz,sr,sv,ampl);
                        //printf("(%d,%d,%d,%d) %f\n", k1,k2,k3,k4,phase[ind(k1,k2,k3,k4)]);
                        //phaseTemp[ind(k1,k2,k3,k4)] = 0;
                            }
                        }
                    }
                }
            }
            //printf("%d\n",k1);
        }

        
    double mass0 = calDensity();
        
    printDensityXY("./datFiles/densXY0.dat",tamano/2);
    printDensityYZ("./datFiles/densYZ0.dat",0);
    printDensityXZ("./datFiles/densXZ0.dat",0);
    
    //printThr(printPotXY,"./datFiles/pot0XY");

    //calAcce();
    //printThr2(printAcceXY,"./datFiles/accex0XY",0);

    potencial(); 

        
  //  printThr(printDensityXY,"./datFiles/densXY");
    //printThr(printPotXY,"./datFiles/pot1XY");
    calAcce();
    
    printThr(printPotXY,"./datFiles/pot0XY");
 //   printThr(printPotXY,"./datFiles/pot1XY");
 //   printThr(printPotXY,"./datFiles/pot1XY");
    
    printThr2(printAcceXY,"./datFiles/accex0XY",0);


    //printAcceXY(char *name, int corteZ, int xyz);
    //printAcceYZ(char *name, int corteX, int xyz);
    //printAcceXZ(char *name, int corteY, int xyz);
    //step();
    //printDensityXY("./datFiles/densXY1.dat",0);
    //printDensityYZ("./datFiles/densYZ1.dat",0);
    //printDensityXZ("./datFiles/densXZ1.dat",0);
    //calDensity();
    
    int suprai = 0;
    printf("La masa es %f\n", mass0);
    
    
    
	for(suprai = 1; suprai<Nt;suprai+=1)
    {    
        char *grid0 = (char*) malloc(200* sizeof(char));
        char *grid1 = (char*) malloc(200* sizeof(char));
        char *grid2 = (char*) malloc(200* sizeof(char));
    
        //printf("Error Mesage00\n");	
		step();
        
        //Descomentar para versión colisional --
       
       // if(TAU != 0){
        //calMacro(); 
        //collisionStep();
        //}
        
        //--
		
        
		sprintf(grid0, "./datFiles/densXY%d.dat", suprai);
        //sprintf(grid1, "./datFiles/densYZ%d.dat", suprai);
        //sprintf(grid2, "./datFiles/densXZ%d.dat", suprai);
        printf("%d %f\n",suprai,calDensity()*100/mass0); //Calcula la densidad.
        printDensityXY(grid0,tamano/2);
        //printDensityYZ(grid1,0);
        //printDensityXZ(grid2,tamano);

        
		potencial();
		//sprintf(grid, "./datFiles/potential%d.dat", suprai);
		//printPot(grid);
        
        
		calAcce();
		//sprintf(grid, "./datFiles/acce%d.dat", suprai);
		//printAcce(grid);
        //sprintf(grid, "./datFiles/gridx%d.dat", suprai);
        //printPhaseX(grid, Ny/2, Nvy/2);
        //sprintf(grid, "./datFiles/gridy%d.dat", suprai);
        //printPhaseY(grid, Nx/2, Nvx/2);
        //printPhase(grid);
        
        free(grid0);
        free(grid1);
        free(grid2);
        
        //step();
        
                
	}

    
    
    
    fclose(constantes);
	//fclose(simInfo);
	return 0;

}

//Calcula el potencial (V) con el método de Fourier. Actualiza el arreglo pot.
double potencial()
{
    inE=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
    inR=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
    mem=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);

    pIda = fftw_plan_dft_3d(Nx, Ny, Nz, inE, out, FFTW_FORWARD, FFTW_MEASURE);

    //Cargar densidad en in:
    for(k1=0;k1<Nx;k1+=1){
        for(k2=0;k2<Ny;k2+=1){
            for(k3=0;k3<Nz;k3+=1){
        inE[in(k1,k2,k3)] = giveDensity(k1,k2,k3) - totalMass/((Xmax-Xmin)*(Ymax-Ymin)*(Zmax-Zmin));
        inR[in(k1,k2,k3)] = 0;
        out[in(k1,k2,k3)] = 0;
            }
        }
    }

    
    
    fftw_execute(pIda);
    
    //Guarda out en mem. 
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2<Ny; k2+=1){
            for(k3 = 0; k3<Nz; k3+=1){
         mem[in(k1,k2,k3)] = out[in(k1,k2,k3)];   
            }
        }
    }

    pIda = fftw_plan_dft_3d(Nx, Ny, Nz, out, inR,FFTW_BACKWARD, FFTW_MEASURE);
    //Se debe usar el mismo plan sí o sí al parecer.

    //Devuelve carga a out Î(Chi).
    
    
    
    //out[0] = -4*PI*G*mem[0] ;
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny   ;k2 += 1){
                for(k3 = 0; k3 <Nz   ;k3 += 1){
        
             out[in(k1,k2,k3)] = -4.0*PI*1.0*(mem[in(k1,k2,k3)])*calcK2((double)k1,(double)k2,(double)k3);//Porque dx = dy y Nx = Ny.

            
        //out[in(k1,k2)] = mem[in(k1,k2)]; //Descomentar esta línea para obtener la distribucion original.
                }
        }
    }

    fftw_execute(pIda);


    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1){
            for(k3 = 0; k3 <Ny;k3 += 1){
                pot[in(k1,k2,k3)] = creal(inR[in(k1,k2,k3)]/Nx/Ny/Nz);//+ totalMass/((Xmax-Xmin)*(Ymax-Ymin));
            }
        }
    }

 
 fftw_free(inE);
 fftw_free(out);
 fftw_free(inR);
 fftw_free(mem);
 fftw_free(out2);
 
 return -1;
 
}

//Calcula los coeficientes de Fourier en la aproximación pseudoespectral.
double calcK2(double i2, double j2, double k2)
{
    if( ( (j2 == 0)   )  && ( (i2 == 0)  ) && ( k2 == 0 ) ){
        return 0;
    }
    if ( (j2 == Ny/2+1)  && (i2 == Nx/2+1) && (k2 == Nz/2+1) ) {
        return 0;
    }

     if(i2<Nx/2){
         i2 = PI*i2;
     }
     if(j2<Nx/2){
         j2 = PI*j2;
     }
     if(k2<Nx/2){
         k2 = PI*k2;
     }
     
     if(i2>=Nx/2){
         i2 = -PI*(Nx-i2);
     }
     if(j2>=Nx/2){
         j2 = -PI*(Ny-j2);
     }    
     if(k2>=Nx/2){
         k2 = -PI*(Nz-k2);
     }    
    return 1.0/(i2*i2+j2*j2+k2*k2);

}

//Imprime las constantes de la simulación para referencia fuera del C.
void printConstant(char *name, double value)
{
    fprintf(constantes, "%s", name);
    fprintf(constantes, " %f\n", value);
}

//Retorna el valor de la gaussiana para un x,v, sigma x, sigma v, y una amplitud dada.
double gaussD(double x,double y, double z, double vx, double vy, double vz, double sr, double sv, double amplitude)
{
	//double ex = -x*x/(sr)-y*y/(sr)-vx*vx/(sv)-vy*vy/(sv);
    	double ex = -x*x/(2.0*sr*sr)-y*y/(2.0*sr*sr)-z*z/(2.0*sr*sr)-vx*vx/(2.0*sv*sv)-vy*vy/(2.0*sv*sv)-vz*vz/(2.0*sv*sv);
        return amplitude*exp(ex);

}

//Calcula la densidad. Actualiza el arreglo density
double calDensity()
{
    totalMass = 0;
        for(k1=0; k1<Nx;k1+=1){
            for(k2=0; k2< Ny; k2+=1){
                for(k3=0; k3< Nz; k3+=1){
                    density[in(k1,k2,k3)] = 0;
                    for(k4=0; k4< Nvx; k4+=1){
                        for(k5=0; k5<Nvy;k5+=1){
                            for(k6=0; k6<Nvz;k6+=1){
                                density[in(k1,k2,k3)] += phase[ind(k1,k2,k3,k4,k5,k6)]*dvx*dvy*dvz;
                                //density[in(k1,k2,k3)] = giveDensity(k1,k2,k3)+ phase[ind(k1,k2,k3,k4,k5,k6)];
                            }
                        }
                    }
                    totalMass += density[in(k1,k2,k3)]*dz*dy*dx;
                }
            }
        }
        return totalMass;
}

//Calcula las variables macroscópicas para la función de equilibrio.
void calMacro()
{
    totalMass = 0;
        for(k1=0; k1<Nx;k1+=1){
            for(k2=0; k2< Ny; k2+=1){
                for(k3=0; k3< Nz; k3+=1){
                    density[in(k1,k2,k3)] = 0;
                    velocityx[in(k1,k2,k3)] = 0;
                    velocityy[in(k1,k2,k3)] = 0;
                    velocityz[in(k1,k2,k3)] = 0;
                    energy[in(k1,k2,k3)] = 0;
                    for(k4=0; k4< Nvx; k4+=1){
                        for(k5=0; k5<Nvy;k5+=1){
                            for(k6=0; k6<Nvz;k6+=1){
                                density[in(k1,k2,k3)] += phase[ind(k1,k2,k3,k4,k5,k6)]*dvx*dvy*dvz;
                                velocityx[in(k1,k2,k3)] += phase[ind(k1,k2,k3,k4,k5,k6)]*dvx*dvy*dvz*darVx(k4);
                                velocityy[in(k1,k2,k3)] += phase[ind(k1,k2,k3,k4,k5,k6)]*dvx*dvy*dvz*darVy(k5);
                                velocityz[in(k1,k2,k3)] += phase[ind(k1,k2,k3,k4,k5,k6)]*dvx*dvy*dvz*darVz(k6);
                                energy[in(k1,k2,k3)] += phase[ind(k1,k2,k3,k4,k5,k6)]*dvz*dvy*dvx*(pow(darVx(k4) - velocityx[in(k1,k2,k3)],2) + pow(darVy(k5) - velocityy[in(k1,k2,k3)],2) + pow(darVz(k6) - velocityz[in(k1,k2,k3)],2))/2.0;
                            }
                        }
                    }
                    totalMass += density[in(k1,k2,k3)]*dz*dy*dx;
                }
            }
        }
}

//Retorna la densidad en (in1,in2,in3).
double giveDensity(int in1, int in2, int in3)
{
    double rta = density[in(in1,in2,in3)];
    return rta;
}

//Retorna el potencial en (in1,in2,in3).
double givePot(int in1, int in2, int in3)
{
    double rta = pot[in(in1,in2,in3)];
    return rta;
}

//Retorna aceleración en (x=0,y=1,z=2) del punto (in1,in2,in3).
double giveAcce(int in1, int in2, int in3, int xyz)
{
    double rta = acce[ina(in1,in2,in3,xyz)];
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
    return -1;
}

//Deriva el potencial y carga la aceleración en el arreglo acce.
void calAcce()
{
    for(k1 = 0; k1<Nx ; k1 +=1){
        for(k2 = 0; k2<Ny ; k2 +=1){
            for(k3 = 0; k3<Nz ; k3 +=1){
    acce[ina(k1,k2,k3,0)] =  (-givePot(mod(k1-2,Nx),k2,k3) + 8*givePot(mod(k1-1,Nx),k2,k3)-8*givePot(mod(k1+1,Nx),k2,k3)+givePot(mod(k1+2,Nx),k2,k3))/(12*dx);
    acce[ina(k1,k2,k3,1)] =  (-givePot(k1,mod(k2-2,Ny),k3) + 8*givePot(k1,mod(k2-1,Ny),k3)-8*givePot(k1,mod(k2+1,Ny),k3)+givePot(k1,mod(k2+2,Ny),k3))/(12*dy);
    acce[ina(k1,k2,k3,2)] =  (-givePot(k1,k2,mod(k3-2,Nz)) + 8*givePot(k1,k2,mod(k3-1,Nz))-8*givePot(k1,k2,mod(k3+1,Nz))+givePot(k1,k2,mod(k3+2,Nz)))/(12*dz);
            }
        }
    }
}

//Calcula la nueva posición de (x,y,vx,vy)=(iinx,iiny,jinx,jiny).
double newij(int iinx, int iiny, int iinz, int jinx, int jiny, int jinz)
{
        
        //double x = Xmin*1.0+dx*iinx; //Inicialización
        deltavx = giveAcce(iinx,iiny, iinz, 0)*dt; //Cambio de velocidad.
        djx = deltavx/dvx;                
        djx = (int)djx;                     //Cambio de vel en el tablero.
        j2x = jinx+djx;                     //Nuevo j.

        //double y = Ymin*1.0+dy*iiny; //Inicialización
        deltavy = giveAcce(iinx,iiny, iinz, 1)*dt;
        djy = deltavy/dvy;
        djy = (int)djy;
        j2y = jiny+djy;

        deltavz = giveAcce(iinx,iiny, iinz, 2)*dt;
        djz = deltavz/dvz;
        djz = (int)djz;
        j2z = jinz+djz;
        
        if(j2x < 0 || j2x >= Nvx || j2y < 0 || j2y >= Nvy || j2z < 0 || j2z >= Nvy)
        {
            totalPerdido += phase[ind(iinx,iiny,iinz,jinx,jiny,jinz)];
            return -1;
        }
        deltavx = Vxmin*1.0+dvx*j2x;
        deltavy = Vymin*1.0+dvy*j2y;
        deltavz = Vzmin*1.0+dvz*j2z;
        
        deltax = deltavx*dt;
        dix = deltax/dx*scale;
        dix = (int) dix;
        
        deltay = deltavy*dt;
        diy = deltay/dy*scale;
        diy = (int) diy;
        
        deltaz = deltavz*dt;
        diz = deltaz/dz*scale;
        diz = (int) diz;
        
        i2x = iinx + dix;
        i2x = mod(i2x,Nx);
        
        i2y = iiny + diy;
        i2y = mod(i2y,Ny);
        
        i2z = iinz + diz;
        i2z = mod(i2z,Nz);
//	printf("%d\n",j2);
    return 0;
}


//Calcula un paso de streaming. Guarda una copia del phase actual en phaseOld. Actualiza phase. 
void step()
{
    for(k1=0;k1<Nx;k1+=1) {
        for(k2=0;k2<Ny;k2+=1) {
            for(k3=0;k3<Nz;k3+=1) {
                for(k4=0;k4<Nvx;k4+=1) {
                    for(k5=0;k5<Nvy;k5+=1) {
                        for(k6=0;k6<Nvz;k6+=1) {
                            if(newij(k1,k2,k3,k4,k5,k6) == 0){
                                phaseTemp[ind(i2x,i2y,i2z,j2x,j2y,j2z)] += phase[ind(k1,k2,k3,k4,k5,k6)];
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    for(k1=0;k1<Nx;k1+=1) {
        for(k2=0;k2<Ny;k2+=1) {
            for(k3=0;k3<Nz;k3+=1) {
                for(k4=0;k4<Nvx;k4+=1) {
                    for(k5=0;k5<Nvy;k5+=1) {
                        for(k6=0;k6<Nvz;k6+=1) {
                                phase[ind(k1,k2,k3,k4,k5,k6)] = phaseTemp[ind(k1,k2,k3,k4,k5,k6)];
                                //printf("phase[%d,%d,%d,%d,%d,%d] = %f]",k1,k2,k3,k4,k5,k6,phase[ind(k1,k2,k3,k4,k5,k6)]);
                                phaseTemp[ind(k1,k2,k3,k4,k5,k6)] = 0;
                            
                        }
                    }
                }
            }
        }
    }
    
} 


//Calcula un paso colisional. Guarda una copia del phase actual en phaseOld. Actualiza phase. 
void collisionStep()
{
    for(k1=0;k1<Nx;k1+=1) {
        for(k2=0;k2<Ny;k2+=1) {
            for(k3=0;k3<Nz;k3+=1) {
                for(k4=0;k4<Nvx;k4+=1) {
                    for(k5=0;k5<Nvy;k5+=1) {
                        for(k6=0;k6<Nvz;k6+=1) {
                            if(newijCol(k1,k2,k3,k4,k5,k6) == 0){
                                phaseTemp[ind(i2x,i2y,i2z,k4,k5,k6)] += collision(k1,k2,k3,k4,k5,k6,TAU) + phase[ind(k1,k2,k3,k4,k5,k6)];
                            }
                        }
                    }
                }
            }
        }
    }
	
    for(k1=0;k1<Nx;k1+=1) {
        for(k2=0;k2<Ny;k2+=1) {
            for(k3=0;k3<Nz;k3+=1) {
                for(k4=0;k4<Nvx;k4+=1) {
                    for(k5=0;k5<Nvy;k5+=1) {
                        for(k6=0;k6<Nvz;k6+=1) {
                                phase[ind(k1,k2,k3,k4,k5,k6)] = phaseTemp[ind(k1,k2,k3,k4,k5,k6)];
                                phaseTemp[ind(k1,k2,k3,k4,k5,k6)] = 0;
                            
                        }
                    }
                }
            }
        }
    }
	
}


//Calcula el cambio en r. Ignora j.
double newijCol(int iinx, int iiny, int iinz, int jinx, int jiny, int jinz)
{        
        j2x = jinx; 
        j2y = jiny; 
        j2z = jinz; 

        if((j2x < 0 || j2x >= Nvx) && (j2y < 0 || j2y >= Nvy) && (j2z < 0 || j2z >= Nvz) ) return -1;

        vx = darVx(j2x);
        vy = darVy(j2y);
        vz = darVz(j2z);
        
        x = vx*dt*scale; 
        x = x/dx;
        x = (int) x;
        
        y = vy*dt*scale; 
        y = y/dy;
        y = (int) y;
        
        z = vz*dt*scale; 
        z = z/dz;
        z = (int) z;
        

        i2x = iinx + x;
        i2x = mod(i2x,Nx);
        
        i2y = iiny + y;
        i2y = mod(i2x,Nx);
        
        i2z = iinz + z;
        i2z = mod(i2z,Nz);
//	printf("%d\n",j2);
    return 0;
}

//Calcula la contribución del término colisional.
double collision(int icolx, int icoly, int icolz,  int jcolx, int jcoly, int jcolz, double tau)
{
    if(TAU==0) return 0;
    return (feq(icolx,icoly,icolz,jcolx,jcoly,jcolz) - phase[ind(icolx,icoly,icolz,jcolx,jcoly,jcolz)])/tau; 
}

//Distribución de equilibrio usada para calcular la contribución del término colisional.
double feq(int iposx, int iposy, int iposz, int jvelx, int jvely, int jvelz)
{
    double ex = -1.0*(pow(darVx(jvelx)-velocityx[in(iposx,iposy,iposz)],2) + pow(darVy(jvely)-velocityy[in(iposx,iposy,iposz)],2) + pow(darVy(jvelz)-velocityy[in(iposx,iposy,iposz)],2) )/(2.0*energy[in(iposx,iposy,iposz)]);
    double other = giveDensity(iposx,iposy,iposz) / pow(2.0*PI*energy[in(iposx,iposy,iposz)],3/2);
    return other * exp(ex);    
}

//Calcula la posición del elemento (in1,in2,in3,in4,in5,in6) del espacio de fase (x,y,z,vx,vy,vz).
int ind(int in1, int in2, int in3, int in4, int in5, int in6)
{

    return in6 + Nvz*(in5 + Nvy*(in4+Nvx*(in3 + Nz*(in2 +Ny*in1))));
}

//Calcula la posición del elemento (in1, in2, in3) en los arreglos bidimensionales.
int in(int in1, int in2, int in3)
{
    return in3 + Nz*(in2 +Ny*in1);
}
//Calcula la posición del elemento (in1, in2, in3) en la dimensión in4 x=0,y=1,z=2. (para cantidades vectorales)
int ina(int in1, int in2, int in3, int in4)
{
    return in4+3*(in3 + Nz*(in2 +Ny*in1));
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

//Métodos para pasar de enteros a valor coordenado.
double darX(int input)
{
 return Xmin*1.0+dx*input;   
}
double darY(int input)
{
 return Ymin*1.0+dy*input;   
}
double darZ(int input)
{
 return Zmin*1.0+dz*input;   
}
double darVx(int input)
{
 return Vxmin*1.0+dvx*input;   
}
double darVy(int input)
{
 return Vymin*1.0+dvy*input;   
}
double darVz(int input)
{
 return Vzmin*1.0+dvz*input;   
}

//Métodos para imprimir los cortes.

//Imprime el plano XY en el corte z=corteZ.
void printDensityXY(char *name, int corteZ)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){
            fprintf(output,"%f ", giveDensity(i,j,corteZ));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}

//Imprime el plano XZ en el corte y=corteY.
void printDensityXZ(char *name, int corteY)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nz;j+=1){ 
            fprintf(output,"%f ", giveDensity(i,corteY,j));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}
//Imprime el plano YZ en el corte x=corteX.
void printDensityYZ(char *name, int corteX)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Ny;i+=1) {
		for(j=0;j<Nz;j+=1){ 
            fprintf(output,"%f ", giveDensity(corteX,i,j));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}

//Imprime el plano XY en el corte z=corteZ de la aceleración xyz. 
void printAcceXY(char *name, int corteZ, int xyz)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
            fprintf(output,"%f ", giveAcce(i,j,corteZ,xyz));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}

//Imprime el plano XY en el corte z=corteZ de la aceleración xyz. 
void printAcceXY2(char *name, int corteZ, int xyz)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
            fprintf(output,"%f ", acceTheo[ina(i,j,corteZ,xyz)]);
        }
		fprintf(output,"\n");
			}
	fclose(output);
}

//Imprime el plano XZ en el corte y=corteY de la aceleración xyz. 
void printAcceXZ(char *name, int corteY, int xyz)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nz;j+=1){ 
            fprintf(output,"%f ", giveAcce(i,corteY,j,xyz));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}
//Imprime el plano YZ en el corte x=corteX de la aceleración xyz. 
void printAcceYZ(char *name, int corteX, int xyz)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Ny;i+=1) {
		for(j=0;j<Nz;j+=1){ 
            fprintf(output,"%f ", giveAcce(corteX,i,j,xyz));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}



//Imprime el plano XY en el corte z=corteZ.
void printPotXY(char *name, int corteZ)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
            fprintf(output,"%f ", givePot(i,j,corteZ));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}

//Imprime el plano XZ en el corte y=corteY.
void printPotXZ(char *name, int corteY)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nz;j+=1){ 
            fprintf(output,"%f ", givePot(i,corteY,j));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}
//Imprime el plano YZ en el corte x=corteX.
void printPotYZ(char *name, int corteX)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Ny;i+=1) {
		for(j=0;j<Nz;j+=1){ 
            fprintf(output,"%f ", givePot(corteX,i,j));
        }
		fprintf(output,"\n");
			}
	fclose(output);
}

//Imprime todos los cortes de una función.
void printThr(void (*f)(char*, int), char *name)
{
    char *copia = (char*) malloc(200* sizeof(char));
    char *num = (char*) malloc(200* sizeof(char));
    for(k1=0; k1<tamano; k1+=1){
    copia = (char*) malloc(200* sizeof(char));
    num = (char*) malloc(200* sizeof(char));
    strcpy(copia,name);
    sprintf(num, "%d.dat", k1);
    strcat(copia,num);
        f(copia,k1);
    }
}

//Imprime todos los cortes de una función en la coordenada dada.
void printThr2(void (*f)(char*, int, int), char *name, int xyz)
{
    char *copia = (char*) malloc(200* sizeof(char));
    char *num = (char*) malloc(200* sizeof(char));
    for(k1=0; k1<tamano; k1+=1){
    copia = (char*) malloc(200* sizeof(char));
    num = (char*) malloc(200* sizeof(char));
    strcpy(copia,name);
    sprintf(num, "%d.dat", k1);
    strcat(copia,num);
        f(copia,k1,xyz);
    }
}





//Funciones de prueba para el poisson solver. TODO: borrar cuando sea seguro.

double darAcceTheo(double x, double y, double z, int xyz)
{
        double rta = 2.0*gaussD2(x,y,z,sr,sv,ampl)/(sv*sv);
        if(xyz == 0)
        {
        //printf("fue x\n");
         return rta * x;
        }
        if(xyz == 1)
        {
         return rta * y;
        }
        if(xyz == 2)
        {
         return rta * z;
        }
    return 0.0;
    
}


double gaussD2(double x,double y, double z, double sr, double sv, double amplitude)
{
    	double ex = -x*x/(sr*sr)-y*y/(sr*sr)-z*z/(sr*sr);
        return amplitude*exp(ex);
}

double gaussD2Dens(double x,double y, double z, double sr, double sv, double amplitude)
{
	//double ex = -x*x/(sr)-y*y/(sr)-vx*vx/(sv)-vy*vy/(sv);
        return gaussD2(x,y,z,sr,sv,amplitude)*(4.0*x*x+4.0*y*y-6.0*sr*sr+4.0*z*z)/(4.0*PI*G*pow(sr,4));

}

double densidadTeorica2(int inx, int iny, int inz, int nx, int ny, int nz)
{
    
 return -(pow(PI*0.5*2.0*nx,2)*sin(PI*darX(inx)*0.5*nx*2.0)+pow(PI*0.5*nz*2.0,2)*sin(PI*darZ(inz)*0.5*nz*2.0)+pow(PI*0.5*ny*2.0,2)*sin(PI*darY(iny)*0.5*ny*2.0))/(4.0*PI);

}

//Potencial con sin(nx*pi*x/Nx) + sin(ny*pi*y/Ny)
double potencialTeorico2(int inx, int iny, int inz, int nx, int ny, int nz)
{
 return sin(PI*darX(inx)*0.5*nx*2)+sin(PI*darY(iny)*0.5*ny*2.0)+sin(PI*darZ(inz)*0.5*nz*2.0);
 //np.sin(np.pi*xin*nx*2.0/L)+np.sin(np.pi*yin*ny*2.0/L)

}









