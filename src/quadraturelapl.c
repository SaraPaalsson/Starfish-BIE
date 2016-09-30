#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>


/*
Computes the solution u to Laplace's equation at all z, using complex density mu. With standard 16-point composite Gauss-Legendre quadrature on panels. 

Note: No z on boundary!
*/
void normalquadlapl(double *u, double *mu, double complex *z, double complex *zDrops, double complex *zpDrops, double *wDrops, int Nz, int Nboundary) {

	unsigned int idom, jbound;

	// Ugly. Fix this.
	const double pi = 3.1415926535897932385;

	for (idom = 0; idom<Nz; idom++) {
		u[idom] = 0*u[idom];
		for (jbound = 0; jbound < Nboundary; jbound++) {
			u[idom] = u[idom] +  mu[jbound]*wDrops[jbound]*cimag(zpDrops[jbound]/(zDrops[jbound]-z[idom]));
		}
		u[idom] = u[idom]*0.5/pi;
	}

}