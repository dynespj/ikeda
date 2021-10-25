#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pari/pari.h>
#include "common.h"





int main (int argc, char **argv)
{
    long i, j; long count;

    /* parse command line */
    if (argc < 4) {
        printf("usage: make-diagonal-matrix <n> <p> <filename> <lambda_1> ... <lambda_n> <filename>\n");
        exit(-1);
    }

    pari_init(100000000,2);
    long n = atol(argv[1]);
    long p = atol(argv[2]);

    if (argc != 4+n) {
        printf("usage: make-diagonal-matrix <n> <p> <filename> <lambda_1> ... <lambda_n> <filename>\n");
        exit(-1);
    }

    long *lambda = (long *) malloc(n*sizeof(long));
    for (i = 0; i < n; i++)
        lambda[i] = atol(argv[4+i]);

    long len = 0; for (i = 1; i <= n; i++) len += i;
    GEN *B = (GEN *) malloc(len*sizeof(GEN));
    GEN gen_p = stoi(p);

    long *ind = (long *) malloc(n*n*sizeof(long));
    count = 0;

    for (i = 0; i < n; i++)
        for (j = i; j < n; j++) {ind[n*i+j] = ind[n*j+i] = count++;}

    for (i = 0; i < len; i++)
        B[i] = gen_0;

    for (i = 0; i < n; i++)
        B[ind[n*i+i]] = gpowgs(gen_p,lambda[i]);

    long *data, datacount;
    pari_pack(B, &data, &datacount, len);

    /* save data to file */
    FILE *fp = fopen(argv[3], "wb");
    long dim = 2; long lens[2] = {n,n};
    fwrite(&dim, sizeof(long), 1, fp);
    fwrite(lens, sizeof(long), 2, fp);
    fwrite(data, sizeof(long), datacount, fp);

    fclose(fp); pari_close();
    return 0;

    return 0;
}
