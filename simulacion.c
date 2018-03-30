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

//Valores límites para la posición.
#define Xmin -1
#define Xmax 1

//Valores límites para la velocidad
#define Vmin -1
#define Vmax 1

//Tamaño del espacio
#define Nx 2048
#define Nv 2048

//Constantes de unidades
#define alpha 5.402
#define aMetros 24
#define aSegundos 18
#define aMasasSol 11


//Arreglos
double phase[Nx][Nv];
double phaseOld[Nx][Nv];
double phaseTemp[Nx][Nv];
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

double dx = (Xmax-Xmin)*1.0/Nx;
double dv = (Vmax-Vmin)*1.0/Nv;
double dt = 0.05; //Se toma 0.5 para repetir los resultados de Franco. 0.5 en mis unidades equivale a ~3mil millones de años. Hay que repensar dispersion de vel.
int Nt = 5;
FILE *constantes;

//Métodos
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


int main()
{
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
	double x;
	double v;
	for(i=0;i<Nx;i+=1) {
		for(j=0;j<Nv;j+=1){
			x = Xmin*1.0+dx*i;
			v = Vmin*1.0+dv*j;
			phase[i][j] = gaussD(x,v,0.1,0.1,4); //1 de dispersion de velocidad equivale a 1000 km/s. Valor tomado del Coma Cluster.
			phaseOld[i][j] = 0;
			phaseTemp[i][j] = 0;
				}
			}
	printPhase("grid1.dat");
	double mass = convertir(calDensity(),aMasasSol)/pow(10,14); // lo divido entre 10**14 para comparar con el coma cluster. En el cual se basó el sistema de unidades.
	printf("%f\n",mass);
	printDensity("density.dat");

    potencial();

    calAcce();
    printAcce("acce.dat");
    step();
	for(int suprai = 0; suprai<Nt;suprai+=1){
		calDensity();
		potencial();
		calAcce();
		step();
	}
    printPhase("grid2.dat");


	fclose(constantes);
	return 0;

}



//Imprime el espacio de fase con el String name como nombre.
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
	return amplitude*exp(ex)/(sx*sv);

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

    pIda = fftw_plan_dft_1d(Nx, out, inR, FFTW_BACKWARD, FFTW_MEASURE); //Se debe usar el mismo plan sí o sí al parecer.
    double memDo;

    //pVuelta2 = fftw_plan_dft_c2r_1d(Nx, out2, inR, FFTW_MEASURE);


    //Devuelve carga a out Î(Chi).
    //out2[0] = mem[0];
    out[0] = mem[0];
    for(i=1;i<Nx;i+=1){
      out[i] = -mem[i]/calcK2((double)i);
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

//Método para evitar efectos misticos relacionado a pointers.
double giveDensity(int l)
{
    double rta = density[l];
    return rta;
}


double convertir(double valor, int unidad )
{
    if(unidad == aMasasSol){
        return valor/(1.988*4*PI*alpha*(6.674)*pow(10,-aMasasSol));
    }
    if(unidad == aMetros){
        return valor/(alpha * pow(10,-aMetros));
    }
    if( unidad == aSegundos){
        return valor/(alpha*pow(10,-aSegundos));
    }
}

void calAcce()
{
    for(i = 0; i<Nx ; i +=1){
    acce[i] =  (pot[(i+1) % Nx] - pot[i])/dx;
    }
}

void printAcce(char *name)
{
	FILE *output = fopen(name, "w+");
	for(i=0;i<Nx;i+=1) {
            fprintf(output, "%f\n",acce[i]);
			}
	fclose(output);
}


//Calcula la posición del bloque i,j en el instante actual+dt. Retorna -1 si se sale del tablero
double newij(int iin, int jin)
{
        double x = Xmin*1.0+dx*iin; //Conversión estandar usada en el main durante la inicialización del phase.
        double v = Vmin*1.0+dv*jin;


        v += acce[iin]*dt;
        x += v*dt;



        j2 = (v-Vmin*1.0)/dx;
        if(j2 <0 || j2 > Nv) return -1;
        i2 = (x-Xmin*1.0)/dv;

        i2 = mod((int) i2,Nx);

	j2 = (int) j2;
//	printf("%d\n",j2);
    return 1;
}

//Calcula un paso. Guarda una copia del phase actual en phaseOld. Actualiza phase. k,i son x. j,l son v.
void step()
{
	for(k = 0; k<Nx; k++){
		for(l= 0; l<Nv; l++){
			if(newij(k,l) >0){
				phaseOld[k][l] = phase[k][l];
				phaseTemp[i2][j2] += phase[k][l];
			}
		}
	}

	for(i = 0; i<Nx; i++){
		for(j= 0; j<Nv; j++){
			phase[i][j] = phaseTemp[i][j];
		}
	}



}

//Observando que el m = p % q es negativo si p<0 y q>0, se define una función de módulo.
int mod(int p, int q)
{
	p = p%q;
	if(p<0){
		 return p+q;
	}
	return p;

}

















