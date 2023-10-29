/* ---------- Atomic Vibrations in 1-D Lattices -----------
Written by Jacksen Narvaez, 2023.
Based on the original code from F. Dominguez, (2000)
---------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N       1000    // Number of atoms in the chain
#define NP      100     // Number of averages
#define NL      1000    // number of frequencies
#define NM      2       // Types of atoms

#define RAN() ((double)rand()/(double)(RAND_MAX))

int concent(int const nm) {
    double dx = 1./(double) nm;
    double rd = RAN();
    return (int) floor(rd/dx);
}

void init_mass_rand (double alpha[], int const n, double const mass[], int const nm) {
    int ii;
    for (ii = 0; ii<n; ii++) {
        alpha[ii] = mass[concent(nm)]/mass[0];
    }
}

void init_mass_per (double alpha[], int const n, double const mass[], int const nm) {
    int ii;
    int jj;
    for (ii = 0; ii<n; ii+=nm) {
        for (jj = 0; jj<nm; jj++) {
            alpha[ii+jj] = mass[jj]/mass[0];
        }
    }
}

int ratio_n(double const alpha[], double const lambda, int const n) {
    double ra = 0.;
    double rd;
    int k, count = 0;
    for(k=0; k<n; k++) {
        rd = 2. - alpha[k]*lambda - ra;
        if (fabs(rd)<1e-50) {
            count += 1;
            k++;
            ra = 0.;
        } else {
            if(rd<0.) count += 1;
            ra = 1./rd;
        }
    }
    return count;
}

int main() {
    double  mass[NM]     = {1., 2.};    // Mass of A atoms
    double  lmin         = 0.;          // Minimum frequency 
    double  lmax         = 5.;          // Maximun frequency 
    char    filename[]   = "./results/Datos2R.txt";
    int     rd           = 1;           // Random or Periodic chain
    double  alpha[N];                   // Reduced masses
    double  idos[NL];                   // Accumulative density modes

    srand(1234);

    int ii, jj;
    double lambda, dlambda = (double)(lmax-lmin)/NL;
    FILE *out;
    out = fopen(filename, "w");
    for(jj=0; jj<NL; jj++) {
        idos[jj] = 0.;
    }
    if (rd == 0) {
    	    init_mass_per(alpha, N, mass, NM);
        }
    for(ii=0; ii<NP; ii++) { 
        if (rd==1) {
    	    init_mass_rand(alpha, N, mass, NM);
        }
        for(jj=0; jj<NL; jj++) {
            lambda = lmin+dlambda*(double)jj;
            idos[jj] += ratio_n(alpha, lambda, N);
        }
    }
    for(jj=0; jj<NL; jj++) {
        fprintf(out, "%8.61f    ", lmin+dlambda*(double)jj);
        fprintf(out, "%8.61f\n", idos[jj]/(double)N/(double)NP);
    }
    fclose(out);

    return 0;
}