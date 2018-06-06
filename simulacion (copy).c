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
#define Vmin -1.0
#define Vmax 1.0

//Tamaño del espacio
#define Nx 2048
#define Nv 2048

//Constantes de unidades - Aún no hay unidades que simulen el Coma Cluster.
#define aMetros 18
#define aSegundos 14
#define aByear 4
#define aMasasSol 5
//#define G 8.5295e-06
#define G 0.000213238135515

#define solarMases 1e14
#define fracT0 5e-2
#define mParsecs 10


//Arreglos
double phase[Nx][Nv] = {0};

double phaseOld[Nx][Nv] = {0};
double phaseTemp[Nx][Nv] = {0};
double *density;
double *pot;
double *acce;

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

double dt = 0.5; //Se toma 0.5 para repetir los resultados de Franco. 0.5 en mis unidades equivale a ~3mil millones de años. Hay que repensar dispersion de vel.
int Nt = 25;
FILE *constantes;
void printPhase(char *name);
double gaussD(double x, double v, double sx, double sv, double amplitude);
double calDensity();
void printDensity(char *name);
void printConstant(char *name, double value);
double giveDensity(int l);
double potencial();
double calcK2(double j2);
double convertir(double valor, int unidad);
void calAcce();
void printAcce(char *name);
double newij(int iin, int jin);
void step();
int mod(int p, int q);
void printPot(char *name);


int main()
{
    dt = dt*dx/dv;
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
	double vSx = 0.1;
	double vSv = 0.05;
	double ampl = 1;
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nv;j+=1){
			x = Xmin*1.0+dx*i;
			v = Vmin*1.0+dv*j;
			phase[i][j] = gaussD(x,v,vSx,vSv,ampl); //1 de dispersion de velocidad equivale a 1000 km/s. Valor tomado del Coma Cluster.
			phaseOld[i][j] = 0;
			phaseTemp[i][j] = 0;
				}
			}
	//printPhase("grid1.dat");
	double mass = calDensity();

	printf("Se simuló %f millones de años con %d pasos de %f millones de años cada uno\n", convertir(Nt,aByear),Nt, convertir(dt,aByear));
	printf("La masa fue %f masas coma\n",convertir(mass, aMasasSol)/1e14);

	printDensity("density.dat");

	FILE *simInfo = fopen("./images/simInfo.dat","w+");
	fprintf(simInfo,"Para la simulación se utilizó las siguientes condiciones:\n");
	fprintf(simInfo,"x va de %.2f a %.2f , v va de %.2f a %.2f\n", Xmin,Xmax,Vmin,Vmax);
	fprintf(simInfo,"Una distribución gaussiana centrada en 0 para el espacio de fase con (sx sv A)=\n");
	fprintf(simInfo,"(%.3f %.3f %.3f)\n", vSx, vSv, ampl);
	fprintf(simInfo,"Se simuló %d instantes con dt = %.3f\n", Nt,dt);

    potencial();

    calAcce();
    printAcce("acce.dat");
    printf("G es %lf\n", G*1.0);


	for(int suprai = 0; suprai<Nt;suprai+=1){
        char *grid = (char*) malloc(200* sizeof(char));
		sprintf(grid, "./datFiles/grid%d.dat", suprai);

        //printf("Error Mesage00\n");
		printPhase(grid);
		step();

		//calDensity();
		printf("%d %f\n",suprai,calDensity()*100/mass);
		sprintf(grid, "./datFiles/density%d.dat", suprai);
		printDensity(grid);

		potencial();
		sprintf(grid, "./datFiles/potential%d.dat", suprai);
		printPot(grid);
		calAcce();
		sprintf(grid, "./datFiles/acce%d.dat", suprai);
		printAcce(grid);
		free(grid);

	}
//    printPhase("grid2.dat");


	fclose(constantes);
	fclose(simInfo);
	return 0;

}



//Imprime el espacio de fase con el String name como nombre.
void printPhase(char *name)
{
	FILE *output = fopen(name, "w+");

	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nv;j+=1){
          //      printf("ignorarPrimero\n");
			fprintf(output,"%f ", phase[i][j]);
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
            fprintf(output, "%f\n",density[i]);
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

//Calcula la densidad. Actualiza el arreglo density
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

//Calcula el potencial (V) con el método de Fourier. Actualiza el arreglo pot.
double potencial()
{
//    double * densityIN= malloc(sizeof(double)*Nx);
//    for(i = 0;i<Nx;i+=1){
//        densityIN[i] = giveDensity(i);
//    }


    fftw_complex *in, *out, *inR, *mem, *out2;
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    inR=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    mem=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);


    //fftw_complex *out2;
   // out2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
    //double *in2;
    //in2 = (double*) malloc((sizeof(double)*Nx));

    fftw_plan pIda, pVuelta;
    pIda = fftw_plan_dft_1d(Nx, in, out,FFTW_FORWARD, FFTW_MEASURE);

    //fftw_plan pIda2, pVuelta2;
    //pIda2 = fftw_plan_dft_r2c_1d(Nx, in2, out2, FFTW_MEASURE);

    FILE *input = fopen("inF.dat", "w+");
    FILE *output0 = fopen("outF0.dat", "w+");
    FILE *output1 = fopen("outF1.dat", "w+");
    FILE *oR = fopen("oR.dat", "w+");
    FILE *oI = fopen("oI.dat", "w+");


    //Cargar densidad en in:
    for(i=0;i<Nx;i+=1){
        //in[0][i] = densityIN[i];
        in[i] = giveDensity(i);
        inR[i] = -1.0;
      //  in2[i] = giveDensity(i);
        //in[i] = sin(2.0*PI*i*deltax);//Actualmente funciona para sin(x) confirmado. (se esperaba que la parte real de out fuera 0 y la imaginaria tuviera los picos, sin embargo solo el módulo cumple esto)
        //printf("%f, %d,%d\n", in[0][i], i,Nx);
        fprintf(input, "%f\n",creal(in[i]));
    }
//    printf("Press enter to continue...\n");
//    getchar();


    fftw_execute(pIda);
//    fftw_execute(pIda2);

    //Guarda out en mem. Imprime a archivo.
    for(i=0;i<Nx;i+=1){
        mem[i] = out[i];
        fprintf(output0, "%f\n",creal(out[i]));
        fprintf(output1, "%f\n",cimag(out[i]));


//        mem[i] = out2[i];
//        fprintf(output0, "%f\n",creal(out2[i]));
//        fprintf(output1, "%f\n",cimag(out2[i]));

    }

    fftw_execute(pIda);
    pIda = fftw_plan_dft_1d(Nx, out, inR, FFTW_BACKWARD, FFTW_MEASURE); //Se debe usar el mismo plan sí o sí al parecer.
    double memDo;

    //pVuelta2 = fftw_plan_dft_c2r_1d(Nx, out2, inR, FFTW_MEASURE);


    //Devuelve carga a out Î(Chi).
    //out2[0] = mem[0];
    out[0] = -4*PI*G*mem[0];
    for(i=1;i<Nx;i+=1){
      out[i] = -4*PI*G*mem[i]/calcK2((double)i);
    //out[i] = mem[i]; //Descomentar esta línea para obtener la distribucion original.
    //out2[i] = mem[i];


      //  memDo = out[i];//Evita problemas de casting.
        //printf("%f %f %d \n",memDo, calcK2((double) i), i);
    }



    fftw_execute(pIda);
    //fftw_execute(pVuelta2);

//        for(i=0;i<Nx;i+=1){
//        fprintf(oR, "%f\n",inR[i]/Nx);
//        fprintf(oI, "%f\n",cimag(inR[i])/Nx);
//    }

    for(i=0;i<Nx;i+=1){
        pot[i] = creal(inR[i]/Nx);
        fprintf(oR, "%f\n",creal(inR[i])/Nx);
        fprintf(oI, "%f\n",cimag(inR[i])/Nx);
    }

    fclose(input);
    fclose(output0);
    fclose(output1);
    fclose(oR);
    fclose(oI);

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

//Calcula la aceleración a partir del arreglo pot actual.

//void calAcce()
//{
//    for(i=1;i<Nx-1;i++){
//        acce[i]=(-pot[i]+pot[i-1])/(dx);
//  }
//    acce[0]=-(pot[0]-pot[Nx-1])/(dx);
//    acce[Nx-1]=-(pot[Nx-1]-pot[Nx-2])/(dx);
//}

//void calAcce()
//{
//    for(i = 0; i<Nx ; i    +=1){
//    acce[i] =  (-pot[mod(i+1,Nx)] + pot[mod(i-1,Nx)])/(2*dx);
//    }
//}

//Deriva el potencial y carga la aceleración en el arreglo acce.
void calAcce()
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
            fprintf(output, "%f\n",acce[i]);
			}
	fclose(output);
}

//Imprime el arreglo Pot
void printPot(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",pot[i]);
			}
	fclose(output);
}



//Version que discretiza los cambios y no el nuevo valor.
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
        x = v*dt;
        double di = x/dx;
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
				phaseOld[k][l] = phase[k][l];
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

//Observando que el m = p % q es negativo si p<0 y q>0, se define una función de módulo con rango de 0 a q-1.
int mod(int p, int q)
{
	p = p%q;
	if(p<0){
		 return p+q;
	}
	return p;

}

















