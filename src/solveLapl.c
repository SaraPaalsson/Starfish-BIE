#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "math.h"

void laplreaddata(int* sizevec, int* Npanels, double complex* panels, double complex *zDrops, double complex *zpDrops, double complex *zppDrops, double* wDrops,double *tpar, double complex *zDom, double *mu, double *unorm, double *uspec, double *ucorrect, char* fileData); 

void normalquadlapl(double *u, double *mu, double complex *z, double complex *zDrops, double complex *zpDrops, double *wDrops, int Nz, int Nboundary);

void specialquadlapl(double *u_specq, double *u_standardq, double *mu, double complex *zDom, double complex *zDrops, double complex *zpDrops, double *wDrops, double complex *panels, int N, int Npanels);

/* 
Compute solution to Laplace's equation with normal and special quadrature given discretization data and complex density as given in laplData.txt.
*/

int main(void){

	FILE *ptr_file;
	char buf[100000];
	int param, rows, tmp;
	int *sizevec,*Npanels;
	unsigned int i,j;
    double complex *zDrops, *zpDrops, *zppDrops, *zDom, *panels;
    double *wDrops, *tpar, *mu, *unorm, *uspec, *ucorrect;
    double *u_standardq, *u_specq;

    double *testzRE, *testzIM, aRE, aIM, bRE, bIM, *errorvec;
    double complex *testz, a, b, tmpAB;
    double umax, errormax;
    char fileData[] = "laplData.txt";
    char *fileDataPtr; 
    fileDataPtr = fileData; 
/*
Read in sizes of variables
*/
	ptr_file = fopen(fileData,"r");
    param = 17;
    rows = 2;
//    int sizevec[param];
    sizevec = malloc(param*sizeof(int));
	if (!ptr_file)
		return 1;
	for (i=0; i < rows; i++) {
		if (i == 1) {
			for (j = 0; j < param; j++) {
				fscanf(ptr_file, "%d",&sizevec[j]);
			}
		} else {
			fgets(buf,100000,ptr_file);
		}
	}
	fclose(ptr_file);

/*
Declare variables for domain, discretization and precomputed values
*/
	Npanels = malloc(sizevec[0]*sizeof(int));
	panels = malloc(sizevec[1]*sizeof(complex double));
	zDrops = malloc(sizevec[3]*sizeof(complex double));
	zpDrops = malloc(sizevec[5]*sizeof(complex double));
	zppDrops = malloc(sizevec[7]*sizeof(complex double));
	wDrops = malloc(sizevec[9]*sizeof(double));
	tpar = malloc(sizevec[10]*sizeof(double));
	zDom = malloc(sizevec[11]*sizeof(complex double));
	mu = malloc(sizevec[13]*sizeof(double));
	unorm = malloc(sizevec[14]*sizeof(double));
	uspec = malloc(sizevec[15]*sizeof(double));
	ucorrect = malloc(sizevec[16]*sizeof(double));

/*
Read in variables from laplData.txt
*/
	laplreaddata(sizevec, Npanels, panels, zDrops, zpDrops, zppDrops, wDrops, tpar, zDom, mu, unorm, uspec, ucorrect, fileDataPtr);


/*
Compute solution with standard quadrature
*/
	u_standardq = malloc(sizevec[11]*sizeof(double));
	normalquadlapl(u_standardq, mu, zDom, zDrops, zpDrops, wDrops, sizevec[11], sizevec[3]);


/*
Correct computed solution using special quadrature where needed.
*/

	u_specq = malloc(sizevec[14]*sizeof(double));
	specialquadlapl(u_specq,u_standardq,mu,zDom,zDrops,zpDrops,wDrops,panels,sizevec[11],sizevec[1]);

//	printf("uspec[i]=%f\n",u_specq[0]);

	umax = 0;
	for (i=0; i<sizevec[14]; i++) {
		if (fabs(umax) < fabs(ucorrect[i])) {
			umax = fabs(ucorrect[i]);
		}
	}
	printf("umax = %f\n",umax );

	errorvec = malloc(sizevec[11]*sizeof(double));
	errormax = 0;
	for (i=0; i<sizevec[14]; i++) {
		errorvec[i] = fabs(u_specq[i]-ucorrect[i])/umax;
		//printf("errorvec[%d]=%f\n",i,errorvec[i] );
		if (errormax < errorvec[i]) {
			errormax = errorvec[i];
		}
	}
	printf("Max error: %12.16f \n", errormax);


/*
	for (i=131; i<132; i++) {
		printf("Obtained unorm[%d]=%12.16f\n",i,unorm[i]);
		printf("Computed u[%d]=%12.16f\n",i,u_standardq[i]);
		printf("Uknown[%d] = %12.16f\n",i, ucorrect[i]);
		printf("uspec[%d] =%12.16f\n", i, uspec[i]);
		printf("Computed u_specq[%d] = %12.16f\n",i,u_specq[i] );
		printf("error[%d] = %12.8f\n",i, errorvec[i] );
	}
*/



	return 0;
}