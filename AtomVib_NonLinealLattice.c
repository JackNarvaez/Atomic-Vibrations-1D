#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N       1000    // Number of atoms in the chain
#define MA      1.0     // Mass of A atoms
#define MB      1.3     // Mass of B atoms
#define C       0.5     // concentration of B atoms
#define NP      100     // Number of averages
#define NL      1000    // number of frequencies
#define LMIN    0.0     // Minimum frequency 
#define LMAX    5.0     // Maximun frequency 
#define FILENAME "DatosNL.txt"

#define RAN() ((double)rand()/(double)(RAND_MAX))

double alphard[N];      // Reduced masses
double idosrd[NL];      // Accumulative density modes

double alphaic[N];      // Reduced masses
double idosic[NL];      // Accumulative density modes

void massrand (void){
    int i;
    for (i = 0; i<N; i++){
        if(RAN()>= C) alphard[i] = 1.;
        else alphard[i] = MB/MA;
    }
}

void massic (void){
    int i;
    for (i = 0; i<N; i+=4){
        alphaic[i] = 1.;
	alphaic[i+1] = MB/MA;
	alphaic[i+2] = MC/MA;
	alphaic[i+1] = MD/MA;
    }
}

int main(void) {
    int i,j,k;
    double rard, raic, rdrd, rdic, dlambda=(double)(LMAX-LMIN)/NL;
    FILE *out;
    time_t t;
    printf("%d  \r", 1);
    srand ((unsigned)time(&t));
    out = fopen(FILENAME, "w");
    for(j=0; j<NL; j++) {
        idosrd[j] = 0.;
        idosic[j] = 0.;
    }
    massic();
    for(i=0; i<NP; i++) { 
    	massrand();       
        for(j=0; j<NL; j++) {
            rard = 0.;
            raic = 0.;
            for(k=0; k<N; k++) {
                rdrd = 2. - alphard[k]*(LMIN+dlambda*(double)j)-rard;
                if (fabs(rdrd)<1e-50) {
                    idosrd[j] += 1;
                    k++;
                    rard = 0.;
                } else {
                    if(rdrd<0.) idosrd[j] += 1.;
                    rard = 1./rdrd;
                }
            }
            for(k=0; k<N; k++) {
                rdic = 2. - alphaic[k]*(LMIN+dlambda*(double)j)-raic;
                if (fabs(rdic)<1e-50) {
                    idosic[j] += 1;
                    k++;
                    raic = 0.;
                } else {
                    if(rdic<0.) idosic[j] += 1.;
                    raic = 1./rdic;
                }
            }
        }
    }
    for(j=0; j<NL; j++) {
        fprintf(out, "%8.61f    ", LMIN+dlambda*(double)j);
        fprintf(out, "%8.61f    ", idosic[j]/(double)N/(double)NP);
        fprintf(out, "%8.61f\n", idosrd[j]/(double)N/(double)NP);
    }
    fclose(out);

    return 0;
}
