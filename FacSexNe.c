/* FacSexNe.c 

Simulation of tracking neutral mutation under facultative sex
To calculate Ne under different scenarios
Code based on that used in Agrawal and Hartfield 2016 (https://github.com/MattHartfield/BalSelSims).

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.

Run by executing:
./FacSexNe N sex gc reps
Where:
- N is the population size
- sex is the frequency of sex (a value between 0 = obligate asex, and 1 = obligate sex)
- gc is the frequency of mitotic gene conversion
- reps is how many times to introduce linked neutral allele

*/

/* Preprocessor statements */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Function prototypes */
void neutinit(double *geninit, unsigned int N);
void reproduction(double *geninit, double sex);
void gconv(double *geninit, double gc);
double ncheck(double *geninit);

/* Main program */
int main(int argc, char *argv[]){
	unsigned int N = 0;				/* Population size */
	unsigned int g, i; 				/* Counters. Reps counter, geno counter */
	unsigned int reps = 0;			/* Length of simulation (no. of introductions of neutral site) */
	double sex = 0;					/* Frequency of sexual reproduction */
	double gc = 0;					/* Frequency of gene conversion */
	double Acheck = 0;				/* Frequency of derived neutral allele (A) after each reproduction */
	double Hsum = 0;				/* Summed heterozygosity over transit time of neutral allele */
	char Hout[128];					/* String to hold filename in */
	FILE *ofp_hs = NULL;			/* Pointer for results output */
	
	/* GSL random number definitions */
	const gsl_rng_type * T; 
	gsl_rng * r;
	
	/* This reads in data from command line. */
	if(argc != 5){
		fprintf(stderr,"Invalid number of input values.\n");
		exit(1);
	}
	N = atoi(argv[1]);
	sex = strtod(argv[2],NULL);
	gc = strtod(argv[3],NULL);
	reps = strtod(argv[4],NULL);
	
	/* Arrays definition and memory assignment */
	double *genotype = calloc(3,sizeof(double));				/* Genotype frequencies */
	unsigned int *gensamp = calloc(3,sizeof(unsigned int));		/* New population samples */
	  
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */

	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	/* Setting up neutral genotype */
	neutinit(genotype,N);
	
	/* Reintroducing neutral genotype, resetting hap sum */	
    Acheck = ncheck(genotype);
    Hsum = Acheck*(1.0-Acheck);
    
    /* Introduce and track neutral mutations 'reps' times */
    g = 0;
    while(g < reps){
    	
    	/* Reproduction routine */
    	reproduction(genotype,sex);
    	
    	/* Gene conversion routine */
       	gconv(genotype,gc);
       	
    	/* Sampling based on new frequencies */
       	gsl_ran_multinomial(r,3,N,genotype,gensamp);
       	for(i = 0; i < 3; i++){
       		*(genotype + i) = (*(gensamp + i))/(1.0*N);
       	}
       	
       	/* Checking state of haplotypes: if A fixed reset so can start fresh next time */
		Acheck = ncheck(genotype);
		Hsum += Acheck*(1.0-Acheck);
       	
       	if(Acheck == 0 || Acheck == 1){
       		sprintf(Hout,"/scratch/mhartfield/CC_FSC_Ne/temp_s%.8lf_gc%.8lf.out",sex,gc);
       		if(g != 0){
				ofp_hs = fopen(Hout,"a");
			}else if(g == 0){
				ofp_hs = fopen(Hout,"w");
			}
			fprintf(ofp_hs,"%.10lf\n",Hsum);
			fclose(ofp_hs);
       		g++;
     		 
     		/* Reintroducing neutral genotype, resetting hap sum */     		
	    	neutinit(genotype,N);
			Acheck = ncheck(genotype);
       		Hsum = Acheck*(1-Acheck);
       	}
    
	}	/* End of simulation */
	
	
	/* Freeing memory and wrapping up */
 	gsl_rng_free(r);
 	free(gensamp);
	free(genotype);
	return 0;
}

/* Setting up neutral polymorphism routine */
void neutinit(double *geninit, unsigned int N){
	*(geninit + 0) = 1 - 1/(1.0*N);
	*(geninit + 1) = 1/(1.0*N);
	*(geninit + 2) = 0;
}	/* End of neutral setup routine */

/* Reproduction routine */
void reproduction(double *geninit, double sex){
	/* Fed-in genotype frequencies (for ease of programming) */
	double gaas, gAas, gAAs;
	/* Genotype frequencies after sex */	
	double gaaSX, gAaSX, gAASX;
	/* Genotype frequencies after asex */	
	double gaaAS, gAaAS, gAAAS;
	
	/* Initial definition of genotypes */
	gaas = *(geninit + 0);
	gAas = *(geninit + 1);
	gAAs = *(geninit + 2);
	
	/* Change in frequencies due to sex */
	gaaSX = sex*(gaas + gAas/2.0)*(gaas + gAas/2.0);
	gAaSX = sex*2.0*(gaas + gAas/2.0)*(gAAs + gAas/2.0);
	gAASX = sex*(gAAs + gAas/2.0)*(gAAs + gAas/2.0);
	
	/* Change in frequencies due to asex */
	gaaAS = gaas*(1 - sex);
	gAaAS = gAas*(1 - sex);
	gAAAS = gAAs*(1 - sex);
	
	/* Combining to give overall frequency change following reproduction */
	*(geninit + 0) = gaaAS + gaaSX;
	*(geninit + 1) = gAaAS + gAaSX;
	*(geninit + 2) = gAAAS + gAASX;
		
}	/* End of reproduction routine */

/* Gene conversion routine */
void gconv(double *geninit, double gc){
	
	/* Fed-in genotype frequencies (for ease of programming) */
	double gaar, gAar, gAAr;
	/* Frequencies after gene conversion */
	double gaagc, gAagc, gAAgc;
	
	/* Initial definition of genotypes */
	gaar = *(geninit + 0);
	gAar = *(geninit + 1);
	gAAr = *(geninit + 2);
	
	/* Gene conversion equations */
	gaagc = gaar + gAar*gc/2.0;
	gAagc = gAar*(1.0-gc);
	gAAgc = gAAr + gAar*gc/2.0;
	
	/* Output */
	*(geninit + 0) = gaagc;
	*(geninit + 1) = gAagc;
	*(geninit + 2) = gAAgc;
	
}	/* End of gene conversion routine */

/* Has neutral allele fixed or not? Measuring its frequency */
double ncheck(double *geninit){
	/* Fed-in genotype frequencies (for ease of programming) */
	double gAas, gAAs;
	double Atot = 0.0;        /* Total frequency of A */
	
	/* Initial definition of genotypes */
	gAas = *(geninit + 1);
	gAAs = *(geninit + 2);
		
	/* Checking frequency of A */
	Atot = gAAs + gAas/2.0;
	return Atot;
	
}	/* End of neutral check routine */

/* End of program */
