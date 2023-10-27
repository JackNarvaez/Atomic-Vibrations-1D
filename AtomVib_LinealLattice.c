#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N           1000
#define LMAX        5.0
#define LMIN        0.
#define DL          0.1
#define PI          3.141592654
#define FILENAME    "DatosL.txt"

int main(void) {
    double ra, rd, lambda, aux, dlambda=DL;
    int k, idos;
    FILE *out;
    out = fopen(FILENAME, "w");
    for (lambda=LMIN; lambda <= LMAX; lambda += dlambda) {
        ra = 0.;
        idos = 0;
        for (k=1; k<N; k++){
            rd = 2.0 -lambda - ra;
            if (fabs(rd)<1e-50) {
                idos ++;
                k++;
                ra = 0.;
            } else {
                if (rd<0.) idos ++;
                ra = 1./rd;
            }
        }
        if (lambda<4.) aux = (2./PI)*asin(sqrt(lambda/4.));
        else aux=1.;
        fprintf(out, "%8.61f    %8.61f  %8.61f \n", lambda, (double)idos/N, aux);
    }
    fclose(out);
    return 0;
}
