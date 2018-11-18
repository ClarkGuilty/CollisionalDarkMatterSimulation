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




#define mParsecs 50e-3  //Cuántos megaparsecs equivalen a una unidad espacial.
#define solarMases 1e12 //Cuántas masas solares equivalen a una unidad de masa.
#define fracT0 3e-3     //Qué fracción de la edad del universo equivale a una unidad de tiempo
#define G 0.006141 //G en estas unidades. Se calcula con calculations.py

#define scale 1.0 //1.0 es el valor estandar, modificarlo para mejorar visualización.

//Unidades funcionales para clusters galácticos.
//#define mParsecs 5
//#define solarMases 1e15
//#define fracT0 2e-1
//#define G 0.2729448134597113


#define TAU 11963


//Arreglos
//static double phase[Nx][Ny][Nvx][Nvy] = {0};
double *phase;
double *phaseTemp;
double *density;
double *velocityx;
double *velocityy;
double *energy;
double *pot;
double *accex;
double *accey;
double *accext;
double *acceyt;

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
int suprai;

double totalMass = 0;

fftw_complex *inE, *out, *inR, *mem, *out2;
fftw_plan pIda;
fftw_plan pVuelta;

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
double dx = (Xmax - Xmin)*1.0/Nx;
double dy = (Ymax - Ymin)*1.0/Nx;
double dvx = (Vxmax - Vxmin)*1.0/Nvx;
double dvy = (Vymax - Vymin)*1.0/Nvy;

double dt = 0.4;
int Nt = 50;

double totalPerdido;


FILE *constantes;

double gaussD(double x,double y, double vx, double vy, double sr, double sv, double amplitude); 
double calDensity();
void calAcce(); 
double potencial(); 
double giveDensity(int in1,int in2);
double giveAccex(int in1, int in2);
double giveAccey(int in1, int in2);
double calcK2(double i2, double j2);
double convertir(double valor, int unidad);
double newij(int iinx, int jinx, int iiny, int jiny); 
void step(); 
int mod(int p, int q);
void printPot(char *name);
int ind(int in1, int in2, int in3, int in4);
int in(int in1, int in2);
double calcK4(double i2, double j2);
double feq2(int ipos, int jvel); //TODO
double feq(int iposx, int iposy, int jvelx, int jvely);
double collision(int icolx, int icoly, int jcolx, int jcoly, double tau);
void collisionStep();
double newijCol(int iinx, int iiny, int jinx, int jiny);
void calMacro();

//Métodos para pasar de enteros a valor coordenado.
double darX(int input);
double darY(int input);
double darVx(int input);
double darVy(int input);

//Métodos de impresión.
void printPhaseX(char *name, int corteY, int corteVy);
void printPhaseY(char *name, int corteX, int corteVx);
void printDensity(char *name);
void printConstant(char *name, double value);
void printAcce(char *namex, char *namey);



int main()
{
    dt = dt*dx*dy/dvx/dvy;
    phase = malloc((sizeof(double)*Nx*Ny*Nvx*Nvy));
    phaseTemp = malloc((sizeof(double)*Nx*Ny*Nvx*Nvy));
    printf("Alloqué memoria exitosamente\n");
    if(phase == NULL){
        printf("phase es Null\n");   
        }
    density = malloc((sizeof(double)*Nx*Ny));
    accex = malloc((sizeof(double)*Nx*Ny));
    accey = malloc((sizeof(double)*Nx*Ny));
    velocityx = malloc((sizeof(double)*Nx*Ny));
    velocityy = malloc((sizeof(double)*Nx*Ny));
    pot = malloc((sizeof(double)*Nx*Ny));
    energy = malloc((sizeof(double)*Nx*Ny));

    
    accext = malloc((sizeof(double)*Nx*Ny*Nt));
    acceyt = malloc((sizeof(double)*Nx*Ny*Nt));
    
	constantes = fopen("./datFiles/constants.dat","w+");
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
    printConstant("TAU", TAU);
	double x;
	double vx;
    double y;
	double vy;
	double sr = 0.2;
    double sv = 0.1;
	double ampl = 5.0;
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
                        vy = Vymin*1.0 + dvy*k4;
                        phase[ind(k1,k2,k3,k4)] = gaussD(x,y,vx,vy,sr,sv,ampl); 
                    }
                }
            }
            //printf("%d\n",k1);
        }
        double mass0 = calDensity();
        printf("Masa = %f\n", mass0);
        //printf("en 64 = %f\n", phase[ind(64,64,64,64)] );
        printDensity("./datFiles/density0.dat");
        printPhaseX("./datFiles/gridx0.dat", Ny/2, Nvy/2);
        printPhaseY("./datFiles/gridy0.dat", Nx/2, Nvx/2);
        potencial();
        printPot("./datFiles/potential0.dat");
        //printf("aun sirve\n");
        calAcce();
        printAcce("./datFiles/accex0.dat", "./datFiles/accey0.dat");
        totalPerdido = 0;
        //step();
        //printPhaseX("./datFiles/gridx1.dat", Ny/2, Nvy/2);
        //printPhaseY("./datFiles/gridy1.dat", Nx/2, Nvx/2);
        //double mass1 = calDensity();
        //printDensity("./datFiles/density1.dat");
        //printf("Masa2 = %f\n", mass1/mass0);
        //printf("Total = %f\n", (mass1+totalPerdido)/mass0);
        //printDensity("Density1.dat");
        printf("Se simuló %f millones de años con %d pasos de %f millones de años cada uno\n", convertir(Nt*dt,aByear)*1000,Nt, convertir(dt,aByear)*1000);
        //fclose(constantes);


        
           printf("G es %lf\n", G*1.0);


	for(suprai = 1; suprai<Nt;suprai+=1){
        char *grid = (char*) malloc(200* sizeof(char));
        //printf("Error Mesage00\n");
		
		step();
        
        //Descomentar para versión colisional --
        if(TAU != 0){
        calMacro(); 
        collisionStep();
        }
        //--
		
        
		sprintf(grid, "./datFiles/density%d.dat", suprai);
        printf("%d %f\n",suprai,calDensity()*100/mass0); //Calcula la densidad.
		printDensity(grid);

		potencial();
		sprintf(grid, "./datFiles/potential%d.dat", suprai);
//		printPot(grid);
		calAcce();
		sprintf(grid, "./datFiles/acce%d.dat", suprai);
		//printAcce(grid);
        sprintf(grid, "./datFiles/gridx%d.dat", suprai);
        printPhaseX(grid, Ny/2, Nvy/2);
        sprintf(grid, "./datFiles/gridy%d.dat", suprai);
        printPhaseY(grid, Nx/2, Nvx/2);
        //printPhase(grid);
        free(grid);
        
        
        //step();
        
                
	}

	
    fclose(constantes);
	//fclose(simInfo);
	return 0;

}



//Imprime el corte y = corteY , Vy = corteVy del espacio de fase (un plano 2d).
void printPhaseX(char *name, int corteY, int corteVy)
{
 
   	FILE *output = fopen(name, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Nvx+1;j+=1){ 
            fprintf(output,"%f ",phase[ind(i,corteY,Nvx-j,corteVy)]);
        }
		fprintf(output,"\n");
			}
	fclose(output);

    
}

//Imprime el corte x = corteX , Vx = corteVx del espacio de fase (un plano 2d).
void printPhaseY(char *name, int corteX, int corteVx)
{
 
   	FILE *output = fopen(name, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=1;j<Nvx+1;j+=1){ 
            fprintf(output,"%f ",phase[ind(corteX,i,corteVx,Nvx-j)]);
        }
		fprintf(output,"\n");
			}
	fclose(output);

    
}
//Imprime el arreglo density con el String name como nombre.
void printDensity(char *name)
{
    FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
            fprintf(output,"%f ", giveDensity(i,j));  //Al parecer esta SÍ es la real.
        }
		fprintf(output,"\n");
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

//Actualiza los arreglos de las variables macroscópicas (densidad, velocidad del fluido y energía interna).
void calMacro()
{
	for(k1=0; k1<Nx;k1+=1){
        for(k2=0; k2< Ny; k2+=1){
            density[in(k1,k2)] = 0;
            velocityx[in(k1,k2)] = 0;
            velocityy[in(k1,k2)] = 0;
            energy[in(k1,k2)] = 0;
            for(k3=0; k3< Nvx; k3+=1){
                for(k4=0; k4< Nvy; k4+=1){
                    density[in(k1,k2)] = giveDensity(k1,k2)+ phase[ind(k1,k2,k3,k4)]*dvy*dvx;
                    velocityx[in(k1,k2)] = velocityx[in(k1,k2)]+ phase[ind(k1,k2,k3,k4)]*dvy*dvx*darVx(k3);
                    velocityy[in(k1,k2)] = velocityy[in(k1,k2)]+ phase[ind(k1,k2,k3,k4)]*dvy*dvx*darVy(k4);
                    energy[in(k1,k2)] = energy[in(k1,k2)]+ phase[ind(k1,k2,k3,k4)]*dvy*dvx*(pow(darVx(k3) - velocityx[in(k1,k2)],2) + pow(darVy(k4) - velocityy[in(k1,k2)],2))/2.0;
                }
            }
            velocityx[in(k1,k2)] = velocityx[in(k1,k2)] / density[in(k1,k2)];
            velocityy[in(k1,k2)] = velocityy[in(k1,k2)] / density[in(k1,k2)];
            energy[in(k1,k2)] = energy[in(k1,k2)] / density[in(k1,k2)];
            totalMass += density[in(k1,k2)]*dx*dy;
        }
    }
    
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
                        density[in(k1,k2)] = giveDensity(k1,k2)+ phase[ind(k1,k2,k3,k4)]*dvy*dvx;
                    }
                }
                totalMass += density[in(k1,k2)]*dx*dy;
            }
        }
        return totalMass;
}

//Calcula el potencial (V) con el método de Fourier. Actualiza el arreglo pot.
double potencial()
{
    inE=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    inR=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    mem=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    pIda = fftw_plan_dft_2d(Nx, Ny, inE, out, FFTW_FORWARD, FFTW_MEASURE);

    //Cargar densidad en in:
    for(k1=0;k1<Nx;k1+=1){
        for(k2=0;k2<Ny;k2+=1){
        inE[in(k1,k2)] = giveDensity(k1,k2) - totalMass/((Xmax-Xmin)*(Ymax-Ymin));
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
    //out[0] = -4*PI*G*mem[0] ;
    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny   ;k2 += 1){
             out[in(k1,k2)] = -4.0*PI*1.0*(mem[in(k1,k2)])*calcK4((double)k1,(double)k2);//Porque dx = dy y Nx = Ny.
        //out[in(k1,k2)] = mem[in(k1,k2)]; //Descomentar esta línea para obtener la distribucion original.
        }
    }
        k1 = 0;
     
    fftw_execute(pIda);
    //printf("Masa total %f\n",totalMass);


    for(k1=0;k1<Nx;k1+=1){
        for(k2 = 0; k2 <Ny;k2 += 1)
        pot[in(k1,k2)] = creal(inR[in(k1,k2)]/Nx/Ny);//+ totalMass/((Xmax-Xmin)*(Ymax-Ymin));

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
    if( ( (j2 == 0) || (j2 == Nx/2) )  && ( (i2 == 0) || (i2 == Nx/2) )  ){
        return 1;
    }
    double rta1= sin(dx*PI*j2);
    double rta2= sin(dx*PI*i2);
    return 1.0/(rta1*rta1+rta2*rta2);
}

double calcK4(double i2, double j2)
{
    if( ( (j2 == 0)   )  && ( (i2 == 0)  )  ){
        return 0;
    }
    if ( (j2 == Nx/2+1)  && (i2 == Nx/2+1)  ) {
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

//Imprime el arreglo Acce.
void printAcce(char *namex, char *namey)
{
	FILE *outputx = fopen(namex, "w+");
    FILE *outputy = fopen(namey, "w+");
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Ny;j+=1){ 
            fprintf(outputx,"%f ", giveAccex(i,j)); 
            fprintf(outputy,"%f ", giveAccey(i,j)); 
        }
		fprintf(outputx,"\n");
        fprintf(outputy,"\n");
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
                    //fprintf(output,"%f ", givePot(i,j));
                    fprintf(output, "%f ",pow(convertir(givePot(i,j), aMetros)/convertir(1.0, aSegundos),2)/givePot(i,j)); //Imprime potencial en J/kg
        }
		fprintf(output,"\n");
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


void collisionStep()
{
    	for(k1 = 0; k1<Nx; k1++){
            for(k2 = 0; k2<Ny; k2++){
                for(k3= 0; k3<Nvx; k3++){
                    for(k4= 0; k4<Nvy; k4++){
                        if(newijCol(k1,k2,k3,k4) ==0){
                            phaseTemp[ind(i2x,i2y,k3,k4)] += collision(k1,k2,k3,k4,TAU) + phase[ind(k1,k2,k3,k4)] ;//+ dt*feq2(k,l)*acce[k]*(darVx(l)-velocity[k])/energy[k];
                        }
                    }
			}
		}
	}

    	for(k1 = 0; k1<Nx; k1++){
            for(k2 = 0; k2<Ny; k2++){
                for(k3= 0; k3<Nvx; k3++){
                    for(k4= 0; k4<Nvy; k4++){
                            phase[ind(k1,k2,k3,k4)] = phaseTemp[ind(k1,k2,k3,k4)] ;
                            phaseTemp[ind(k1,k2,k3,k4)] = 0;
                }
			}
		}
	}
	
}


//Calcula el cambio en r. Ignora j.
double newijCol(int iinx, int iiny, int jinx, int jiny)
{
        double x = Xmin*1.0+dx*iinx; //Inicialización
        double vx = accex[iinx]*dt;
        
        double y = Ymin*1.0+dy*iiny; //Inicialización
        double vy = accey[iiny]*dt;

        
        j2x = jinx; 
        j2y = jiny; 

        if((j2x < 0 || j2x >= Nvx) && (j2y < 0 || j2y >= Nvy) ) return -1;

        vx = darVx(j2x);
        vy = darVy(j2y);
        
        x = vx*dt*scale; 
        double dix = x/dx;
        dix = (int) dix;
        
        y = vy*dt*scale;
        double diy = y/dy;
        diy = (int) diy;
        

        i2x = iinx + dix;
        i2x = mod(i2x,Nx);
        
        i2y = iiny + diy;
        i2y = mod(i2y,Ny);
//	printf("%d\n",j2);
    return 0;
}

//Calcula la contribución colisional en phase[icol][jcol] con un Tau dado.
double collision(int icolx, int icoly, int jcolx, int jcoly, double tau)
{
    if(TAU==0) return 0;
    double df = (feq(icolx,icoly,jcolx,jcoly) - phase[ind(icolx,icoly,jcolx,jcoly)])/tau;    
    return df;
    return 0;
}

double feq(int iposx, int iposy, int jvelx, int jvely)
{
    double ex = -1.0*(pow(darVx(jvelx)-velocityx[in(iposx,iposy)],2) + pow(darVy(jvely)-velocityy[in(iposx,iposy)],2) )/(2.0*energy[in(iposx,iposy)]);
    double other = giveDensity(iposx,iposy) / (2.0*PI*energy[in(iposx,iposy)]);
    return other * exp(ex);    
}

double feq2(int ipos, int jvel)
{
    //double ex = -1.0*pow(darVx(jvel),2)/(2.0*energy[ipos]);
    //double other = density[ipos] / sqrt(2*PI*energy[ipos]);
    //double lowMach = 1.0 + darVx(jvel)*velocity[ipos]/energy[ipos] + pow(darVx(jvel)*velocity[ipos],2)/(2.0*energy[ipos]) - pow(velocity[ipos],2)/(2.0*energy[ipos]);
    //return other * exp(ex)* lowMach;
    return 0;
    
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
double darVx(int input)
{
 return Vxmin*1.0+dvx*input;
}
double darVy(int input)
{
 return Vymin*1.0+dvy*input;   
}






