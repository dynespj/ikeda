#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pari/pari.h>
#include "common.h"





int main (int argc, char **argv)
{
    long i;

    /* parse command line */
    if (argc != 3) {
        printf("usage: diagonalize-over-Zp <p> <filename>\n");
        exit(-1);
    }

    pari_init(100000000,2);
    long p = atol(argv[1]);
    char *fn = argv[2];

    /* read date from file */
    FILE *fp = fopen(fn, "rb");
    struct stat filestat;
    stat(fn, &filestat);

    long ndims;
    fread(&ndims, sizeof(long), 1, fp);
    long *dims = (long *) malloc(ndims*sizeof(long));
    fread(dims, sizeof(long), ndims, fp);
    long n = dims[0];
    

    long datacount = (filestat.st_size-(ndims+1)*sizeof(long))/sizeof(long);
    long *data = (long *) malloc(datacount*sizeof(long));
    fread(data, sizeof(long), datacount, fp);

    /* relink data */
    long len = 0; for (i = 1; i <= n; i++) len += i;

    GEN *A = (GEN *) malloc(len*sizeof(GEN));
    pari_relink(A, data, len);

    // printmatrix(A, n, ind);
    GEN gen_p = stoi(p); GEN D;
    long count = 0;

    pari_sp ltop = avma;
    for (i = 0; i < 20000000; i++) {
        diagonalize_Zp(A, &D, gen_p, n);
        count++;
        if (count > 100000) {
            pari_sp lbot = avma;
            gerepile(ltop, lbot, NULL);
            ltop = avma;
            count = 0;
        }
    }

    // GEN det = gen_1;
    //for (i = 0; i < n; i++) {
        //pari_printf("%Ps\n", D[i]);
        //det = gmul(det, gel(D,i));
    //}
    //pari_printf("\nDet = %Ps\n", det);

    /* clean up */
    free(data); fclose(fp);
    pari_close();
    
    return 0;
}