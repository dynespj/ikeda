#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pari/pari.h>
#include "common.h"





int main (int argc, char **argv)
{
    long i,j;

    /* validate command line arguments */
    if (argc != 2) {
        printf("usage: print-matrix <filename>\n");
        exit(-1);
    }

    /* initialize pari library */
    MPI_Init(NULL, NULL);
    pari_init(100000000,2);

    /* open file for reading */
    FILE *fp = fopen(argv[1], "rb");
    struct stat filestat;
    stat(argv[1], &filestat);

    /* read in dimensions */
    long ndims;
    fread(&ndims, sizeof(long), 1, fp);
    long *dims = (long *) malloc(ndims*sizeof(long));
    fread(dims, sizeof(long), ndims, fp);
    long n = dims[0];
    
    long len = 0; for (i = 1; i <= n; i++) len += i;

    /* read in data block */
    long datacount = (filestat.st_size-(ndims+1)*sizeof(long))/sizeof(long);
    long *data = (long *) malloc(datacount*sizeof(long));
    fread(data, sizeof(long), datacount, fp);

    /* relink data */
    GEN *A = (GEN *) malloc(len*sizeof(GEN));
    pari_relink(A, data, len);

    long *ind = (long *) malloc(n*n*sizeof(long));
    long count = 0;

    for (i = 0; i < n; i++)
        for (j = i; j < n; j++) {ind[n*i+j] = ind[n*j+i] = count++;}

    /* print matrix */
    count = 0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            pari_printf("%Ps \t", A[ind[n*i+j]]);
        pari_printf("\n");
    }

    /* clean up */
    free(data); fclose(fp);
    pari_close();
    MPI_Finalize();

    return 0;
}
