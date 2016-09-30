#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

/* 
Read in text file as input for testing. Input written in MATLAB. 
Contains variables (in order):
	Npanels panels zDropsRE zDropsIM zpDropsRE zpDropsIM 
	     zppDropsRE zppDropsIM wDrops tpar zDomRE zDomIM mu unorm uspec ucorrect  
First line, variable names. 
Second line, length of arrays. 
*/

int main() 
{
	FILE *ptr_file;
    char buf[1000];
    unsigned int i, j; 
    int param, size;
    double *zDropsRE, *zDropsIM, *zpDropsRE, *zpDropsIM, *zppDropsRE, *zppDropsIM;
    double *panels, *wDrops, *tpar, *zDomRE, *zDomIM, *mu, *unorm, *uspec, *ucorrect;
    double *Npanels;
    double complex *zDrops, *zpDrops, *zppDrops, *zDom;


    double tmp; 

    ptr_file = fopen("laplData.txt","r");
    param = 16;
    size = 18;

    int sizevec[param];

	if (!ptr_file)
		return 1;

	for (i = 0; i < size; i++) {
		if (i == 1) {
			for (j = 0; j < param; j++) {
				fscanf(ptr_file, "%d",&sizevec[j]);
			}
		}	
		if (i == 2) {
			Npanels = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &Npanels[j]);
			}
		}
		else if (i == 3) {
			panels = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &panels[j]);
			}
		}
		else if (i == 4) {
			zDropsRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zDropsRE[j]);
			}
		}
		else if (i == 5) {
			zDropsIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zDropsIM[j]);
			}
		}
		else if (i == 6) {
			zpDropsRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zpDropsRE[j]);
			}
		}
		else if (i == 7) {
			zpDropsIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zpDropsIM[j]);
			}
		}
		else if (i == 8) {
			zppDropsRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zppDropsRE[j]);
			}
		}
		else if (i == 9) {
			zppDropsIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zppDropsIM[j]);
			}
		}
		else if (i == 10) {
			wDrops = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &wDrops[j]);
			}
		}
		else if (i == 11) {
			tpar = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &tpar[j]);
			}
		}
		else if (i == 12) {
			zDomRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zDomRE[j]);
			}
		}
		else if (i == 13) {
			zDomIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &zDomIM[j]);
			}
		}
		else if (i == 14) {
			mu = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &mu[j]);
			}
		}
		else if (i == 15) {
			unorm = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &unorm[j]);
			}
		}
		else if (i == 16) {
			uspec = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &uspec[j]);
			}
		}
		else if (i == 17) {
			ucorrect = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf", &ucorrect[j]);
			}
		}
		else {
			fgets(buf,1000,ptr_file);
		}
	}
	

	zDrops = malloc(sizevec[2]*sizeof(complex double));
	for (i=0; i<sizevec[2]; i++) {
		zDrops[i] = zDropsRE[i] + zDropsIM[i]*I;
	}


	for (i=0;i<param;i++) {
		printf("Size is %d\n",sizevec[i]);
	}

	printf("Npanels = %lf\n",Npanels[0]);

	printf("5 first of panels = ");
	for (i=0; i<4; i++) {
		printf("%12.8f", panels[i]);
	}
	printf("\n");


	printf("5 first of zDropsRE = ");
	for (i=0; i<4; i++) {
		printf("%12.8f", zDropsRE[i]);
	}
	printf("\n");

	printf("5 first of zDropsIM = ");
	for (i=0; i<4; i++) {
		printf("%12.8f", zDropsIM[i]);
	}
	printf("\n");

	printf("5 first of ucorrect = ");
	for (i=0; i<4; i++) {
		printf("%12.8f", ucorrect[i]);
	}
	printf("\n");



	fclose(ptr_file);
	return 0;
}
