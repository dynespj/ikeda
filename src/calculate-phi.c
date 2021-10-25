#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pari/pari.h>
#include "common.h"

#define BLOCK_LOW(i,p,n)   ((i)*(n)/(p))
#define BLOCK_HIGH(i,p,n)  (BLOCK_LOW((i)+1,p,n)-1)
#define BLOCK_SIZE(i,p,n)  (BLOCK_LOW((i)+1,p,n)-BLOCK_LOW(i,p,n))
#define BLOCK_OWNER(i,p,n) (((p)*((i)+1)-1)/(n))



bool
next_part(long *part, long dim)
{
    long temp, k = dim-1; 
    while(part[k] == 0) k--;

    if (k > 0) {
        part[k-1]++; 
        temp = part[k]-1;
        part[k] = 0; 
        part[dim-1] = temp;
        return true;
    }

    return false;
}

bool
next_elem(long *elem, long *mod, long dim)
{
    long k = dim-1; 
    while (k > -1)
        if (elem[k] >= mod[k]-1) k--;
        else break;

    if (k > -1) {
        elem[k] += 1; k++;
        while (k < dim) 
        {elem[k] = 0; k++;}
        return true;
    }

    return false;
}

void
calculate_work(long rank, long nprocs, long p, long dim, long nu, long *lambda, long *elem, long *mod, long *seq_size)
{
    long j, k, count; long *pow = (long *) malloc(dim*sizeof(long));

    long *ind = (long *) malloc(dim*dim*sizeof(long)); count = 0;
    for (j = 0; j < dim; j++)
        for (k = j; k < dim; k++) 
            {ind[dim*j+k] = ind[dim*k+j] = count; count++;}

    memset(lambda, 0, dim*sizeof(long)); lambda[dim-1] = nu;
    count = 0;

    do {
        long term = 1;
        for (j = 0; j < dim; j++)
            for (k = 0; k < j*lambda[j]; k++)
                term *= p;
        count += term;
    } while (next_part(lambda,dim));

    long seq_start = BLOCK_LOW(rank,nprocs,count);
    *seq_size = BLOCK_SIZE(rank,nprocs,count);

    memset(lambda, 0, dim*sizeof(long)); lambda[dim-1] = nu;
    count = 0;

    do {
        for (j = 0; j < dim; j++) {
            pow[j] = 1;
            for (k = 0; k < lambda[j]; k++) 
                pow[j] *= p;
        }
        long term = 1;
        for (j = 0; j < dim; j++)
            for (k = j+1; k < dim; k++)
                term *= pow[k];
        if (count+term > seq_start) break;
        count += term;
    } while (next_part(lambda,dim));

    long rem = seq_start-count;
    for (j = 0; j < dim; j++) {  elem[ind[dim*j+j]] = 0;
                                 mod[ind[dim*j+j]] = 1; }

    for (j = 0; j < dim; j++)
        for (k = j+1; k < dim; k++)
            mod[ind[dim*j+k]] = pow[k];

    long slice_size = 1;
    for (j = 0; j < dim; j++)
        for (k = j+1; k < dim; k++)
            slice_size *= pow[k];

    for (j = 0; j < dim; j++) {
        for (k = j+1; k < dim; k++) {
            slice_size /= mod[k];
            long quot = rem/slice_size;
            elem[ind[dim*j+k]] = quot;
            rem = rem - quot*slice_size;
        }
    }

    free(ind); free(pow);
}

void
calculate_Hnp(long p, long dim, GEN **Hnp, long **data, long *maxdeg)
{
    long i, j, k, l; GEN gen_p = stoi(p);
    GEN minpol = gtopoly(mkvecn(3,gen_1,gen_0,gneg(gen_p)),1);
    GEN qf_0 = mkpolmod(gtopoly(mkvecn(1,gen_0),1),minpol);
    GEN qf_1 = mkpolmod(gtopoly(mkvecn(1,gen_1),1),minpol);
    GEN qf_m1 = mkpolmod(gtopoly(mkvecn(1,gen_m1),1),minpol);
    GEN qf_p = mkpolmod(gtopoly(mkvecn(1,gen_p),1),minpol); 
    GEN qf_psqrtp = mkpolmod(gtopoly(mkvecn(2,gen_1,gen_0),1),minpol);
    GEN qf_msqrtp = gneg(qf_psqrtp); GEN qf_pmsqrtp[] = {qf_psqrtp,qf_msqrtp};

    GEN *H = (GEN *) malloc(2*(dim+1)*sizeof(GEN));
    H[0] = H[1] = H[2] = H[3] = H[4] = H[5] = gtopoly(mkvecn(1,qf_1),0);

    for (i = 3; i <= dim; i += 2) {
        GEN quadfact = gtopoly(mkvecn(3,gneg(gpowgs(qf_p,2*((i-1)/2)-1)),qf_0,qf_1),0);
        H[2*i] = gmul(quadfact,H[2*(i-2)]); H[2*i+1] = H[2*i+2] = H[2*i+3] = H[2*i]; 
        }

    for (i = 2; i <= dim; i += 2)
        for (j = 0; j <= 1; j++) {
            GEN linfact = gtopoly(mkvecn(2,gmul(qf_pmsqrtp[j],gpowgs(qf_p,(i-1)/2)),qf_1),0);
            H[2*i+j] = gmul(linfact,H[2*i+j]); }

    GEN polfact = gtopoly(mkvecn(3,qf_m1,qf_0,qf_1),0);
    for (i = 0; i < 2*(dim+1); i++) H[i] = gmul(polfact,H[i]);

        
    *maxdeg = 3+2*((dim-1)/2);
    GEN *V = (GEN *) malloc(4*(dim+1)*(*maxdeg+1)*sizeof(GEN));

    for (i = 0; i <= dim; i++)
        for (j = 0; j <= 1; j++) {
            GEN pol = H[2*i+j]; long deg = degree(pol);
            for (k = 0; k <= *maxdeg; k++) {
                for (l = 0; l <= 1; l++) {
                    long index = 4*i*(*maxdeg+1)+2*(*maxdeg+1)*j+2*k+l; 
                    V[index] = gen_0;
                    if (k <= deg) {
                        GEN term = ((GEN *) truecoeff(pol,k))[2];
                        if (l <= degree(term))
                            V[index] = truecoeff(term,l);
                    }
                }
            }
        }

    long datacount;
    pari_pack(V, data, &datacount, 4*(dim+1)*(*maxdeg+1));
    pari_relink(V, *data, 4*(dim+1)*(*maxdeg+1));
    *Hnp = V; free(H);
}

void
calculate_TGinv(long dim, GEN Tdiag, GEN *G, GEN *TGinv)
{
    long i, j, k, count; long *ind; GEN term, *pow, *B;
    pow = (GEN *) malloc(dim*sizeof(GEN));

    count = 0;
    ind = (long *) pari_malloc(dim*dim*sizeof(long)); count = 0;
    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++) 
            {ind[dim*i+j] = ind[dim*j+i] = count; count++;}

    long len = count;
    B = (GEN *) malloc(len*sizeof(GEN));

    for (i = 0; i < dim; i++)
        B[ind[dim*i+i]] = pow[i] = G[ind[dim*i+i]];

    for (j = 1; j < dim; j++)
        for (i = 0; i < j; i++) {
            B[ind[dim*i+j]] = G[ind[dim*i+j]];
            for (k = i+1; k < j; k++) {
                term = gdiv(gmul(B[ind[dim*i+k]],G[ind[dim*k+j]]),pow[k]);
                B[ind[dim*i+j]] = gsub(B[ind[dim*i+j]],term);
            }
        }

    for(i = 0; i < dim; i++)
        for (j = i; j < dim; j++) {
            TGinv[ind[dim*i+j]] = gen_0;
            for (k = 0; k <= i; k++) {
                term = gmul(B[ind[dim*k+i]],B[ind[dim*k+j]]);
                term = gmul(term,gel(Tdiag,k+1));
                term = gdiv(term,gsqr(pow[k]));
                TGinv[ind[dim*i+j]] = gadd(TGinv[ind[dim*i+j]],term);
            }
            TGinv[ind[dim*i+j]] = gdiv(TGinv[ind[dim*i+j]],gmul(pow[i],pow[j]));
        } 

    free(pow); free(ind);
}

bool
is_p_integral(long dim, GEN *A, GEN p)
{
    long i, j, count; long *ind;

    count = 0;
    ind = (long *) pari_malloc(dim*dim*sizeof(long)); count = 0;
    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++) 
            {ind[dim*i+j] = ind[dim*j+i] = count; count++;}

    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++)
            if (Q_pval(A[ind[dim*i+j]],p) < 0)
                return false;

    return true;
}

void
calculate_p_local_parameters(long dim, GEN *A, GEN p, long *sp, long *lambdap)
{
    GEN Adiag; long i, j, nnz, nz = 0; GEN det = Fp_red(gen_1,p);
    long count, len; long *ind;

    count = 0;
    ind = (long *) pari_malloc(dim*dim*sizeof(long)); count = 0;
    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++) 
            {ind[dim*i+j] = ind[dim*j+i] = count; count++;}

    len = count;
    GEN *Amodp = (GEN *) malloc(len*sizeof(GEN));

    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++) 
            Amodp[ind[dim*i+j]] = Fp_red(A[ind[dim*i+j]],p);

    free(ind);

    diagonalize_Fp(Amodp, &Adiag, p, dim);
    
    for (i = 1; i <= dim; i++)
        if (gequal0(gel(Adiag,i))) nz++;
        else det = gmul(det, gel(Adiag,i));

    nnz = dim - nz; *sp = nz;
    det = gmul(gpowgs(Fp_red(gen_m1,p),nnz/2),det);


    if (nnz > 0) {
        if (nnz % 2 == 1) *sp += (nnz-1)/2;
        else if (kronecker(det,p) == 1) *sp += nnz/2;
        else *sp += nnz/2-1;
    }
   
    if (*sp == dim) *lambdap = 0;
    else *lambdap = 1;

    free(Amodp);
}

void
print_upper_triangular_matrix(GEN *A, long dim, long *ind)
{
    long i, j;

    for (i = 0; i < dim; i++) {
        for (j = 0; j < i; j++) pari_printf("0\t");
        for (j = i; j < dim; j++) pari_printf("%Ps\t", A[ind[i*dim+j]]);
        pari_printf("\n");
    }
    pari_printf("\n");
}

void
print_symmetric_matrix(GEN *A, long dim, long *ind)
{
    long i, j;

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) pari_printf("%Ps\t", A[ind[i*dim+j]]);
        pari_printf("\n");
    }
    pari_printf("\n");
}

int
main (int argc, char **argv)
{
    int nprocs, rank;
    pari_sp ltop, lbot;
    long i, j, k, count, subcount; 

    /* validate command line arguments */
    if (argc != 4) {
        printf("usage: calculate-phi <matrix> <p> <beta>\n");
        exit(-1);
    }

    /* initialize */
    MPI_Init(NULL, NULL);
    pari_init(100000000, 2);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char *fn = argv[1];
    FILE *fp = fopen(argv[1], "rb");
    struct stat filestat;
    stat(fn, &filestat);
    long p = atol(argv[2]);
    long beta = atol(argv[3]);

    long ndims;
    fread(&ndims, sizeof(long), 1, fp);
    long *dims = (long *) malloc(ndims*sizeof(long));
    fread(dims, sizeof(long), ndims, fp);
    long dim = dims[0];
    free(dims);

    long *lambda = (long *) malloc(dim*sizeof(long));
    GEN *pow = (GEN *) malloc(dim*sizeof(GEN));
    long len = 0; for (i = 1; i <= dim; i++) len += i;
    long *elem = (long *) malloc(len*sizeof(long));
    long *mod = (long *) malloc(len*sizeof(long));
    GEN *T = (GEN *) malloc(len*sizeof(GEN));
    GEN *G = (GEN *) malloc(len*sizeof(GEN));
    GEN *TGinv = (GEN *) malloc(len*sizeof(GEN)); 
    GEN Tdiag; GEN gen_p = stoi(p);

    count = 0;
    long *ind = (long *) pari_malloc(dim*dim*sizeof(long)); count = 0;
    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++) 
            {ind[dim*i+j] = ind[dim*j+i] = count; count++;}

    /* read in matrix */
    long datacount = (filestat.st_size-(ndims+1)*sizeof(long))/sizeof(long);
    long *data = (long *) malloc(datacount*sizeof(long));
    fread(data, sizeof(long), datacount, fp);
    pari_relink(T, data, len);

    /* diagonalize T */
    diagonalize_Zp(T, &Tdiag, gen_p, dim);

    /* compute polynomials */
    long *block, maxdeg; GEN *Hnp;
    calculate_Hnp(p, dim, &Hnp, &block, &maxdeg);

    //GEN S[] = {gen_0,gen_0};
    ltop = avma; 

    /* main loop */
    for (i = 0; i <= beta/2; i++) {
        if (beta-2*i <= maxdeg) {
            /* split up work */
            long seq_size;
            calculate_work(rank, nprocs, p, dim, i, lambda, elem, mod, &seq_size);

            if (rank == 0)
                printf("seq_size = %ld\n", seq_size);

            /* inner loop */
            count = 0; subcount = 0; bool done = false;

            do {
                if (count > 0) {
                    memset(elem, 0, len*sizeof(long));

                    for (j = 0; j < dim; j++) 
                        pow[j] = gpowgs(gen_p,lambda[j]);

                    for (j = 0; j < dim; j++)
                        mod[ind[dim*j+j]] = 1;

                    for (j = 0; j < dim; j++)
                        for (k = j+1; k < dim; k++)
                            mod[ind[dim*j+k]] = pow[k];
                }
                //for (j = 0; j < dim; j++) 
                    //pow[j] = gpowgs(gen_p,lambda[j]);

                // reset elem!

                do {
                    if (subcount > 1000) {
                        lbot = avma;
                        gerepile(ltop, lbot, NULL);
                        //gerepileall(ltop, lbot, 2, S, S+1);
                        subcount = 0; ltop = avma;
                    }
                    if (count >= seq_size) {done = true; break;}

                    /* construct G */
                    for (j = 0; j < dim; j++)
                        G[ind[dim*j+j]] = pow[j];

                    for (j = 0; j < dim; j++)
                        for (k = j+1; k < dim; k++)
                            G[ind[dim*j+k]] = stoi(elem[ind[dim*j+k]]);
                    
                    pari_printf("G = \n");
                    print_upper_triangular_matrix(G, dim, ind);

                    /* Calculate T[G^-1] */
                    pari_printf("TGinv = \n");
                    calculate_TGinv(dim, Tdiag, G, TGinv);

                    print_symmetric_matrix(TGinv,dim,ind);

                    /* Add term */
                    if (is_p_integral(dim, TGinv, gen_p)) {
                        long sp, lambdap;
                        calculate_p_local_parameters(dim, TGinv, gen_p, &sp, &lambdap);
                        //for (l = 0; l <= 1; l++)
                            //S[l] = gadd(S[l], Hnp[4*sp*(maxdeg+1)+2*(maxdeg+1)*lambdap+2*(beta-2*i)+l]);
                    }

                    count++; subcount++; ;
                } while(next_elem(elem, mod, len) && !done);
            } while (next_part(lambda,dim) && !done);
        }
    }

    // pari_printf("%Ps %Ps\n", S[0], S[1]);

    /* clean up */
    free(lambda); free(pow); 
    free(elem); free(mod);
    pari_close();
    MPI_Finalize();

    return 0;
}
