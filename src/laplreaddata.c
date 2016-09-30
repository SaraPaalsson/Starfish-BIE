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

void laplreaddata(int* sizevec, int* Npanels, double complex* panels, double complex *zDrops, double complex *zpDrops, double complex *zppDrops, double *wDrops, double *tpar, double complex *zDom, double *mu, double *unorm, double *uspec, double *ucorrect, char* fileData) 
{
	FILE *ptr_file;
    char buf[1000];
    unsigned int i, j; 
    int param, size;
    double *zDropsRE, *zDropsIM, *zpDropsRE, *zpDropsIM, *zppDropsRE, *zppDropsIM, *zDomRE, *zDomIM, *panelsRE, *panelsIM;
    double *Npanelstmp;
    double tmp; 

    ptr_file = fopen(fileData,"r");
    param = 17;
    size = 19;


	if (!ptr_file) {
		printf("Error in file read\n");
		// Add break for error here!
	}
		
	for (i = 0; i < size; i++) {
		if (i == 2) {
			for (j = 0; j < sizevec[i-2]; j++) {
				Npanelstmp = malloc(sizevec[0]*sizeof(double));
				fscanf(ptr_file, "%lf", &Npanelstmp[j]);
			}
		}
		else if (i == 3) {
			panelsRE = malloc(sizevec[1]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &panelsRE[j]);			
			}
		}
		else if (i == 4) {
			panelsIM = malloc(sizevec[2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &panelsIM[j]);
			}
		}
		else if (i == 5) {
			zDropsRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zDropsRE[j]);
			}
		}
		else if (i == 6) {
			zDropsIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zDropsIM[j]);
			}
		}
		else if (i == 7) {
			zpDropsRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zpDropsRE[j]);
			}
		}
		else if (i == 8) {
			zpDropsIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zpDropsIM[j]);
			}
		}
		else if (i == 9) {
			zppDropsRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zppDropsRE[j]);
			}
		}
		else if (i == 10) {
			zppDropsIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zppDropsIM[j]);
			}
		}
		else if (i == 11) {
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &wDrops[j]);
			}
		}
		else if (i == 12) {
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &tpar[j]);
			}
		}
		else if (i == 13) {
			zDomRE = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zDomRE[j]);
			}
		}
		else if (i == 14) {
			zDomIM = malloc(sizevec[i-2]*sizeof(double));
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &zDomIM[j]);
			}
		}
		else if (i == 15) {
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &mu[j]);
			}
		}
		else if (i == 16) {
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &unorm[j]);
			}
		}
		else if (i == 17) {
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &uspec[j]);
			}
		}
		else if (i == 18) {
			for (j = 0; j < sizevec[i-2]; j++) {
				fscanf(ptr_file, "%lf12.16", &ucorrect[j]);
			}
		}
		else {
			fgets(buf,1000,ptr_file);
		}
	}
	

	for (i=0; i<sizevec[0]; i++) {
		Npanels[i] = (int) Npanelstmp[i];
	}

	for (i=0; i<sizevec[1]; i++) {
		panels[i] = panelsRE[i] + panelsIM[i]*I;
	}


	for (i=0; i<sizevec[3]; i++) {
		zDrops[i] = zDropsRE[i] + zDropsIM[i]*I;
		zpDrops[i] = zpDropsRE[i] + zpDropsIM[i]*I;
		zppDrops[i] = zppDropsRE[i] + zppDropsIM[i]*I;
	}

	for (i=0; i<sizevec[11]; i++) {
		zDom[i] = zDomRE[i] + zDomIM[i]*I;
	}
	

/*
	for (i=0;i<param;i++) {
		printf("Size is %d\n",sizevec[i]);
	}

	printf("Npanels = %d\n",Npanels[0]);
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
*/


	fclose(ptr_file);
}
