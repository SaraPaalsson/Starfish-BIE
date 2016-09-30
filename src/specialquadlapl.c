#include <string.h>
#include "Complex.h"
#include <stdio.h>
#include <stdlib.h>


const double pi = 3.1415926535897932385;
const int startsize = 32;
void vandernewton(double complex *T, double complex *b, int N);
void vandernewtonT(double *T, double *b, int N);
void IPmultR(double complex* in, double complex* out);

//16-point Gauss-Legendre weights
static const double W16[16] = {0.0271524594117541,
0.0622535239386479,0.0951585116824928,0.1246289712555339,0.1495959888165767,
0.1691565193950025,0.1826034150449236,0.1894506104550685,0.1894506104550685,
0.1826034150449236,0.1691565193950025,0.1495959888165767,0.1246289712555339,
0.0951585116824928,0.0622535239386479,0.0271524594117541};

//32-point Gauss-Legendre weights
static const double W32[32] = {0.0070186100094701,
0.0162743947309057,0.0253920653092621,0.0342738629130214,0.0428358980222267,
0.0509980592623762,0.0586840934785355,0.0658222227763618,0.0723457941088485,
0.0781938957870703,0.0833119242269467,0.0876520930044038,0.0911738786957639,
0.0938443990808046,0.0956387200792749,0.0965400885147278,0.0965400885147278,
0.0956387200792749,0.0938443990808046,0.0911738786957639,0.0876520930044038,
0.0833119242269467,0.0781938957870703,0.0723457941088485,0.0658222227763618,
0.0586840934785355,0.0509980592623762,0.0428358980222267,0.0342738629130214,
0.0253920653092621,0.0162743947309057,0.0070186100094701};

//16 to 32 point interpolation coefficients. 1st set.
static const double IP1[128] = {0.7082336805923402,
0.4174603456614055,0.1201306210679926,-0.0237384214490220,-0.0247115099296993,
0.0087677652494516,0.0103484922682165,-0.0048829375742948,-0.0058171385997042,
0.0033716100497611,0.0039063895362322,-0.0026884914506444,-0.0029826204379799,
0.0023966133657277,0.0025278996374981,-0.0023520880120266,-0.3533989766039463,
0.1289141029247849,0.4904381439593198,0.4476888554483156,0.1595815187430263,
-0.0431996043883891,-0.0450542761365297,0.0198432780500969,0.0226533071187764,
-0.0127661666173026,-0.0145091849486018,0.0098524148293410,0.0108283966165110,
-0.0086459140673153,-0.0090837161649102,0.0084360245293210,0.2547888140892103,
-0.0794386511305539,-0.1800216753917482,0.1127069946981099,0.4529903213227301,
0.4571030934901456,0.1720456302558366,-0.0572686901163087,-0.0572070079705632,
0.0298642753857761,0.0323480604597623,-0.0212795957958572,-0.0228944350026966,
0.0180248076268620,0.0187756290920826,-0.0173659080015727,-0.1926427644762529,
0.0575395983713891,0.1185502002165091,-0.0605587419142252,-0.1375504134801631,
0.1105111980442395,0.4385063216111399,0.4620135345715323,0.1767505762588966,
-0.0700000771781777,-0.0663627648029001,0.0404551610934199,0.0415343932539965,
-0.0317585071715892,-0.0325161000341283,0.0298340633769956,0.1433033540000506,
-0.0420553811514164,-0.0836233891757775,0.0401361142969071,0.0816535055732694,
-0.0527563287484864,-0.1144994775853014,0.1112128815258087,0.4303680392263895,
0.4662432580208928,0.1793331307963325,-0.0837301037551319,-0.0757157523283246,
0.0540255806818779,0.0532626953367605,-0.0480491590014860,-0.0997192378404317,
0.0290155272553322,0.0567456713515015,-0.0265072815309246,-0.0516841720948336,
0.0312108723382504,0.0603950802150088,-0.0468713962117200,-0.0965418713375640,
0.1113929985433533,0.4229333917674608,0.4726563404781661,0.1839456727063571,
-0.1024390131967896,-0.0906825961622589,0.0783222643612123,0.0589566269168384,
-0.0170808697411231,-0.0331317609494943,0.0152767692550160,0.0292158689453183,
-0.0171477684792897,-0.0317987518605555,0.0230873922781534,0.0424630971897936,
-0.0391477381829754,-0.0774309538666544,0.1069678364314479,0.4098502582609164,
0.4887570918888832,0.2021622184713724,-0.1455583641763743,-0.0195214966778084,
0.0056453278101818,0.0109121889216969,-0.0050042888041768,-0.0094951190796480,
0.0055107724940780,0.0100569812321846,-0.0071340625232678,-0.0126690018860246,
0.0110418399786724,0.0197819310583685,-0.0222335618307413,-0.0445659130687800,
0.0796393408723433,0.3555539698235838,0.5967331669239306};

//16 to 32 point interpolation coefficients. 2nd set
static const double IP2[128] = {0.7138621264850029,
0.4158614649995460,0.1171390533885666,-0.0224309414525152,-0.0223867275212341,
0.0075268332425387,0.0083097853751947,-0.0036134992927350,-0.0038983391485322,
0.0020027759055358,0.0020013610561088,-0.0011449345392067,-0.0010004418238182,
0.0005796227498290,0.0003691229777510,-0.0001148410895266,-0.3731117376310109,
0.1345146978688941,0.5009196712113994,0.4431061809431234,0.1514292542409961,
-0.0388453473755598,-0.0378952348465420,0.0153814075225089,0.0159014848426175,
-0.0079431247886893,-0.0077862576814684,0.0043949156603830,0.0038044673660606,
-0.0021902526754258,-0.0013893468075581,0.0004314370407234,0.2935334077534948,
-0.0904491991507758,-0.2006375430165828,0.1217267282571176,0.4690505693866059,
0.4485149826814497,0.1579049657964081,-0.0484399253991464,-0.0438186361102529,
0.0202761928802727,0.0189425113789293,-0.0103579732562729,-0.0087773455062908,
0.0049826169161143,0.0031336115858597,-0.0009691268935199,-0.2543216125483122,
0.0750746089078695,0.1514059982920834,-0.0749489083470473,-0.1632097247818790,
0.1242574593626827,0.4611915989585189,0.4478105302111243,0.1551400216478456,
-0.0544610911859422,-0.0445314841467965,0.0225651764329818,0.0182471281387628,
-0.0100600543577763,-0.0062187415151346,0.0019078707296151,0.2312943045160101,
-0.0670850646237963,-0.1305709521800638,0.0607297941305177,0.1184505233059086,
-0.0725218317815637,-0.1472268604148448,0.1317870430598982,0.4618288268307383,
0.4434844549025628,0.1471232281426801,-0.0570984665414406,-0.0406678216286402,
0.0209226990236611,0.0124538953296838,-0.0037566466272948,-0.2171239069846040,
0.0624388430107100,0.1195285512638427,-0.0541067920843720,-0.1011439299315847,
0.0578788932057221,0.1047623469874682,-0.0749282556026716,-0.1397580556688505,
0.1429367300147700,0.4680721498192956,0.4348189017231927,0.1332828758133183,
-0.0535184789189454,-0.0286039577327160,0.0082607580054626,0.2087875429056409,
-0.0597829885222319,-0.1135080588890587,0.0507179130294734,0.0929917302083465,
-0.0517207938386748,-0.0897133328585768,0.0600282764461917,0.0999806752076895,
-0.0817026011454862,-0.1393794337995033,0.1600513710164597,0.4830068086060422,
0.4153122181040704,0.1037159232533378,-0.0249698016066503,-0.2049002094488522,
0.0585617629264828,0.1108029670564811,-0.0492413053477593,-0.0895742689267964,
0.0492637410706939,0.0840953326985925,-0.0549762659931624,-0.0884105585949814,
0.0683011464015245,0.1055383029995325,-0.0985990125544629,-0.1556639994547752,
0.2005703021781758,0.5406401699812300,0.3033999036738149};


void specialquadlapl(double *u_specq, double *u_standardq, double *mu, double complex *zDom, double complex *zDrops, double complex *zpDrops, double *wDrops, double complex *panels, int N, int Npanels) {

	int i, k, j;

	double complex nzpan[16], tz[16], tzp[16];
	double complex tz32[32], tzp32[32], nzpan32[32], p32[33], orig32[32];
	double tmpT[16], tmpb[16], tW32[32], tW[32];
	double complex zk, mid, len, nz, lg1, lg2, testsum, o32sum;
	double oldsum, test, plen, newsum, sign;
	int furthercheck;

	double complex *tf32, *tf, new1, modif;
	tf32 = malloc(32*sizeof(complex double));
	tf = malloc(16*sizeof(complex double));

/* Go through all points in domain and compute new solution if needed */
	for (i = 0; i < N; i++ ) {
		u_specq[i] = u_standardq[i];
		zk = zDom[i];

		for (k = 0; k<Npanels-1; k++) {

			mid = (panels[k+1]+panels[k])/2;
			len = (panels[k+1]-panels[k]);

	//Is target point zk close enough to panel k to warrant further tests?
            if (cabs(zk-mid) < cabs(len)) {

            if (i==199) {
                printf("i=199, k=%d \n",k);
            }

            	nz = 2*(zk-mid)/len;
            	oldsum = 0;
            	testsum = 0;
            	lg1 = clog(1-nz);
            	lg2 = clog(-1-nz); 

    //Save panel points and map s.t. real part on interval -1,1
            	for (j=0; j<16; j++) {
            		tz[j] = zDrops[k*16+j];
            		nzpan[j] = 2*(tz[j]-mid)/len;
            	}
	//Is the point between the panel and the real axis?
            	if (creal(nz) > -1 && creal(nz) < 1) {

                   if(cimag(nz) > 0) {
	       			//Above the real axis, check if nz is enclosed by the panel and the real axis.
             		    furthercheck = 0;
            		    for (j=0; j<16; j++) {
            			    if (cimag(nzpan[j])>cimag(nz)) {
            				    furthercheck = 1;
            				    break;
            			    }
            		    }
            		    if (furthercheck) {
            			    for (j=0; j<16; j++) {
            				    tmpT[j] = creal(nzpan[j]);
            				    tmpb[j] = cimag(nzpan[j]);
            			    }
            			    vandernewtonT(tmpT,tmpb,16);
            			    test = tmpb[15];
            			    for (j=14; j>=0; j--) {
            				    test = test*creal(nz) + tmpb[j];
            			    }
            			    if (test > cimag(nz)) {
            				    // Yes, iti s, assuming a reasonably well refined mesh
            				    //Correct the value of the integral
            				    lg1 -= pi*I;
            				    lg2 += pi*I;
            		  	    }
                   		}
                	}
                
            	
            	   if (cimag(nz) < 0 ) {
                    //Below the real axis, check if nz is enclosed by the panel and the real axis
    				    furthercheck = 0;
    				    for (j=0; j<16; j++) {
    					   if (cimag(nzpan[j])<cimag(nz)) {
    						    furthercheck = 1;
    						    break;
    					    }
    				    }        		
    				    if (furthercheck) {
    					   for (j=0; j<16; j++) {
    						    tmpT[j] = creal(nzpan[j]);
    						    tmpb[j] = cimag(nzpan[j]);
    					    }
    					
    					    vandernewtonT(tmpT,tmpb,16);
    					    test = tmpb[15];
    					
    					    for (j=14; j>=0; j--) {
    						    test = test*creal(nz) + tmpb[j];
    					    }
    					
    					    if (test < cimag(nz)) {
    						// Yes, it is, assuming a reasonably well refined mesh
    						// Correct the value of the integral
    						    lg1 += pi*I;
    						    lg2 -= pi*I;
    					    }
    					}
                    }
            	} 

            	p32[0] = lg1-lg2; //Analytical expression for p_0

            	oldsum = 0;
            	testsum = 0;
            	//Compute old contribution to u from panel with standard quadrature
 	           for (j=0; j<16;j++) {
    	        	tzp[j] = zpDrops[k*16+j];
        	    	tf[j] = mu[k*16+j];
            		tW[j] = wDrops[k*16+j];
            		oldsum += tW[j]*tf[j]*cimag(tzp[j]/(tz[j]-zk))/(2*pi);
            		testsum += tW[j]*tzp[j]/(tz[j]-zk);
            	}

    //Does standard quadrature suffice? In that case, we need not do anything.
           		if (cabs(p32[0]-testsum) > 1e-13) {
            	//If not, we first attempt 32-point quadrature.
                //Interpolate boundary and density to 32 point Gauss-Legendre points.

					IPmultR(tf,tf32);
					IPmultR(tz,tz32);
					IPmultR(tzp,tzp32);


					plen = tW[0]/W16[0];
					o32sum = 0;
					for (j=0; j<32; j++) {
						tW32[j] = W32[j]*plen;
						orig32[j] = tW32[j]/(tz32[j]-zk);
						o32sum += tzp32[j]*orig32[j];
					}

					if (cabs(o32sum-p32[0]) < 1e-13) {
						//32 point quadrature suffices

						newsum = 0;
						for (j=0; j<32; j++) {
							newsum += tW32[j]*tf32[j]*cimag(tzp32[j]/(tz32[j]-zk))*0.5/pi;
						}

						u_specq[i] += creal(newsum-oldsum);

					} else {
    	                //Straight up 32 point quadrature doesn't suffice.
        	            //We use interpolatory quadrature instead.

						IPmultR(nzpan,nzpan32);

        				sign = -1;
                    	//Compute the analytic values of the integral of
                    	//1/(tau-z) multiplied by monomials.
                    	for (j=1; j<33; j++) {
                    		p32[j] = nz*p32[j-1]+(1.0-sign)/j;
                    		sign = -sign;
                    	}
                    
	                    //Solve the vandermonde systems to get the quadrature 	weights.
	                    vandernewton(nzpan32,p32,32);

    	                //Compute the correct value of the stokeslet
        	            //integral and subtract the old, inaccurate value.

            	        new1 = 0;
                	    for (j=0; j<32; j++) {
                    		new1 += cimag(p32[j]*tf32[j])*0.5/pi;
                    	}
                    	modif = new1-oldsum;
                    	u_specq[i] += creal(modif); 

					}
           		}
            }
		}
	}
}



void vandernewton(double complex *T, double complex *b, int N) {
	int k, j;

	for (k=1; k<N; k++) {
		for (j=N-1; j>=k; j--) {
			b[j] -= T[k-1]*b[j-1];
		} 
	}

	for (k=N-1; k>=1; k--) {
		for (j=k; j<N; j++) {
			b[j] = b[j]/(T[j]-T[j-k]);
			b[j-1] -= b[j];
		}
	}

}


  /*  void vandernewton(Complex *T, Complex *b, int N) {
        for(int k = 1;k < N;k++)
            for(int j = N-1;j>=k;j--) {
                b[j] -= T[k-1]*b[j-1];
//         b[j+16] -= T[k-1]*b[j+15];
            }
        
        for(int k = N-1;k >= 1;k--)
            for(int j=k;j<N;j++) {
                b[j] /= T[j]-T[j-k];
//         b[j+16] /= T[j]-T[j-k];
                b[j-1] -= b[j];
//         b[j+15] -= b[j+16];
            }
    }*/


void vandernewtonT(double *T, double *b, int N) {
	int k,j;

	for (k=0; k<N-1; k++) {
		for (j=N-1; j>=k+1;j--) {
			b[j] = (b[j]-b[j-1])/(T[j]-T[j-k-1]);
		}
	}

	for (k=N-1; k>=0; k--) {
		for (j=k; j<N-1; j++) {
			b[j] -= T[k]*b[j+1];
		}
	}

}

/*
    void vandernewtonT(double *T, double *b, int N) {
        for(int k = 0;k < N-1;k++)
            for(int j = N-1;j>=k+1;j--)
                b[j] = (b[j]-b[j-1])/(T[j]-T[j-k-1]);
        
        for(int k = N-1;k >= 0;k--)
            for(int j=k;j<N-1;j++)
                b[j] -= T[k]*b[j+1];
    }
*/


void IPmultR(double complex *in, double complex *out) {
	int i, j, ptr;
	double complex t1, t2;

    for(i = 0; i < 16; i++) {
        t1 = 0;
        t2 = 0;
        ptr = i;
        for(j = 0; j < 8; j++) {
            t1 += IP1[ptr]*(in[j]+in[15-j]);
            t2 += IP2[ptr]*(in[j]-in[15-j]);
            ptr += 16;
        }
           
        out[i] = t1+t2;
        out[31-i] = t1-t2;
    }
}

/*
    void IPmultR(Complex* in,Complex* out) {
        
        for(int i = 0;i<16;i++) {
            Complex t1 = 0;
            Complex t2 = 0;
            int ptr = i;
            for(int j = 0;j<8;j++) {
                t1 += IP1[ptr]*(in[j]+in[15-j]);
                t2 += IP2[ptr]*(in[j]-in[15-j]);
                ptr += 16;
            }
            out[i] = t1+t2;
            out[31-i] = t1-t2;
        }
        
    }
    */

