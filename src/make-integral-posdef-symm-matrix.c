#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pari/pari.h>
#include "common.h"





int main (int argc, char **argv)
{
    long i, j, k; long count;

    /* parse command line */
    if (argc != 3) {
        printf("usage: make-integral-posdef-symm-matrix <n> <filename>\n");
        exit(-1);
    }

    pari_init(100000000,2);
    long n = atol(argv[1]);

    /* make triangular matrix */
    srand(time(NULL));


    GEN *A = (GEN *) malloc(n*n*sizeof(GEN));

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            if (j > i) A[n*i+j] = gen_0;
            else A[n*i+j] = stoi(rand()%32-16);
        }
    
    /* make pos-def symm matrix */
    long len = 0; for (i = 1; i <= n; i++) len += i;
    GEN *B = (GEN *) malloc(len*sizeof(GEN));
    count = 0;

    for (i = 0; i < n; i++)
        for (j = 0; j <= i; j++) {
            B[count] = gen_0;
            for (k=0; k < n; k++)
                B[count] =gadd(B[count],gmul(A[n*i+k],A[n*k+j]));
            count++;
        }

    long *ind = (long *) malloc(n*n*sizeof(long));
    count = 0;

    for (i = 0; i < n; i++)
        for (j = 0; j <= i; j++) {ind[n*i+j] = ind[n*j+i] = count++;}

    /* print results */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            pari_printf("%Ps \t", B[ind[n*i+j]]);
        pari_printf("\n");
    }

    /* pack random matrix */
    long *data, datacount;
    pari_pack(B, &data, &datacount, len);

    /* save data to file */
    FILE *fp = fopen(argv[2], "wb");
    long dim = 2; long lens[2] = {n,n};
    fwrite(&dim, sizeof(long), 1, fp);
    fwrite(lens, sizeof(long), 2, fp);
    fwrite(data, sizeof(long), datacount, fp);

    fclose(fp); pari_close();
    return 0;
}
