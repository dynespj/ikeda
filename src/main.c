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
    long temp, i = dim-1;
    while(part[i] == 0) i--;

    if (i > 0) {
        part[i-1]++; 
        temp = part[i]-1;
        part[i] = 0; 
        part[dim-1] = temp;
        return true;
    }

    return false;
}

bool
next_elem(long *elem, long *mod, long len)
{
    long i = len-1; 
    while (i > -1)
        if (elem[i] == mod[i]-1) i--;
        else break;

    if (i > -1) {
        elem[i] += 1; i++;
        while (i < len) 
        {elem[i] = 0; i++;}
        return true;
    }

    return false;
}

void
calculate_Hnp(GEN gen_p, long dim, GEN **Hnp, long **data, long *maxdeg)
{
    long i, j, k;
    pari_sp ltop = avma;

    /* initialize field elements */
    GEN minpol = gtopoly(mkvecn(3,gen_1,gen_0,gneg(gen_p)),1);
    GEN qf_0 = mkpolmod(gtopoly(mkvecn(1,gen_0),1),minpol);
    GEN qf_1 = mkpolmod(gtopoly(mkvecn(1,gen_1),1),minpol);
    GEN qf_m1 = mkpolmod(gtopoly(mkvecn(1,gen_m1),1),minpol);
    GEN qf_p = mkpolmod(gtopoly(mkvecn(1,gen_p),1),minpol); 
    GEN qf_psqrtp = mkpolmod(gtopoly(mkvecn(2,gen_1,gen_0),1),minpol);
    GEN qf_msqrtp = gneg(qf_psqrtp); GEN qf_pmsqrtp[] = {qf_psqrtp,qf_msqrtp};

    /* compute the polynomials */
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

    /* extract the coefficients */
    *maxdeg = 3+2*((dim-1)/2);
    GEN *V = (GEN *) malloc(2*(dim+1)*(*maxdeg+1)*sizeof(GEN));

    for (i = 0; i <= dim; i++)
        for (j = 0; j <= 1; j++) {
            GEN pol = H[2*i+j]; long deg = degree(pol);
            for (k = 0; k <= *maxdeg; k++) {
                long index = 2*i*(*maxdeg+1)+(*maxdeg+1)*j+k;
                V[index] = gen_0;
                if (k <= deg) {
                    GEN term = ((GEN *) truecoeff(pol,k))[2];
                    V[index] = truecoeff(term,k%2);
                }
            }
        }

    /* copy data to heap */
    long datacount;
    pari_pack(V, data, &datacount, 2*(dim+1)*(*maxdeg+1));
    pari_relink(V, *data, 2*(dim+1)*(*maxdeg+1));

    pari_sp lbot = avma;
    gerepile(ltop, lbot, NULL);
    *Hnp = V; free(H);
}

void
calculate_work(long rank, long nprocs, long p, long dim, long nu, long *lambda, long *elem, long *seq_size, long *ind)
{
    long j, k; long count;

    long *pow = (long *) malloc(dim*sizeof(long));

    memset(lambda, 0, dim*sizeof(long)); lambda[dim-1] = nu;
    count = 0;

    /* calculate total work */
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
        count += term;
    } while (next_part(lambda,dim));

    /* calculate individual work */
    long seq_start = BLOCK_LOW(rank,nprocs,count);
    *seq_size = BLOCK_SIZE(rank,nprocs,count);

    /* calculate initial index */
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
    for (j = 0; j < dim; j++) 
        elem[ind[dim*j+j]] = 0;

    long slice_size = 1;
    for (j = 0; j < dim; j++)
        for (k = j+1; k < dim; k++)
            slice_size *= pow[k];

    for (j = 0; j < dim; j++) {
        for (k = j+1; k < dim; k++) {
            slice_size /= pow[k];
            long quot = rem/slice_size;
            elem[ind[dim*j+k]] = quot;
            rem = rem - quot*slice_size;
        }
    }

    free(pow);
}

void
populate_G(GEN *G, long dim, long *pow, long *elem, long *ind)
{
    long i, j;

    for (i = 0; i < dim; i++)
        G[ind[dim*i+i]] = stoi(pow[i]);

    for (i = 0; i < dim; i++)
        for (j = i+1; j < dim; j++)
            G[ind[dim*i+j]] = stoi(elem[ind[dim*i+j]]);
}

void
calculate_TInvG(GEN *G, GEN DiagT, long dim, GEN *TInvG, long *ind)
{
    long i, j, k; GEN term;

    long len = 0;
    for (i = 1; i <= dim; i++) len += i;
    GEN *B = (GEN *) malloc(len*sizeof(GEN));

    for (i = 0; i < dim; i++)
        B[ind[dim*i+i]] = G[ind[dim*i+i]];

    for (j = 1; j < dim; j++)
        for (i = 0; i < j; i++) {
            B[ind[dim*i+j]] = gneg(G[ind[dim*i+j]]);
            for (k = i+1; k < j; k++) {
                term = gmul(B[ind[dim*i+k]],G[ind[dim*k+j]]);
                term = gdiv(term,G[ind[dim*k+k]]);
                B[ind[dim*i+j]] = gsub(B[ind[dim*i+j]],term);
            }
        }

    for(i = 0; i < dim; i++)
        for (j = i; j < dim; j++) {
            TInvG[ind[dim*i+j]] = gen_0;
            for (k = 0; k <= i; k++) {
                term = gmul(B[ind[dim*k+i]],B[ind[dim*k+j]]);
                term = gmul(term,gel(DiagT,k+1));
                term = gdiv(term,gsqr(G[ind[dim*k+k]]));
                TInvG[ind[dim*i+j]] = gadd(TInvG[ind[dim*i+j]],term);
            }
            term = gmul(G[ind[dim*i+i]],G[ind[dim*j+j]]);
            TInvG[ind[dim*i+j]] = gdiv(TInvG[ind[dim*i+j]],term);
        } 
}

bool
is_p_integral(GEN *A, GEN gen_p, long dim, long *ind)
{
    long i, j;

    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++)
            if (Q_pval(A[ind[dim*i+j]], gen_p) < 0)
                return false;

    return true;
}

void
reduce_modp(GEN *TInvG, GEN *ModPTInvG, GEN gen_p, long dim, long *ind)
{
    long i, j;

    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++)
            ModPTInvG[ind[dim*i+j]] = Rg_to_Fp(TInvG[ind[dim*i+j]], gen_p);
}

void
calculate_p_local_parameters(GEN Diag, GEN gen_p, long dim, long *sp, long *lambdap)
{
    long i, nnz, nz = 0; 
    GEN fp_1 = Rg_to_Fp(gen_1,gen_p);
    GEN fp_m1 = Rg_to_Fp(gen_m1,gen_p);
    GEN det = fp_1;
    
    for (i = 1; i <= dim; i++)
        if (gequal0(gel(Diag,i))) nz++;
        else det = gmul(det, gel(Diag,i));

    nnz = dim - nz; *sp = nz;
    det = gmul(gpowgs(fp_m1,nnz/2),det);

    if (nnz > 0) {
        if (nnz % 2 == 1) *sp += (nnz-1)/2;
        else if (kronecker(det,gen_p) == 1) *sp += nnz/2;
        else *sp += nnz/2-1;
    }
   
    if (*sp == dim) *lambdap = 0;
    else *lambdap = 1;
}


int
main (int argc, char **argv)
{
    long i, j, k, l, count, subcount;  bool done;
    if (argc != 4) {printf("usage: calculate-phi <matrix> <p> <beta>\n"); exit(-1);}

    /* initialize */
    MPI_Init(NULL, NULL); pari_init(100000000, 2);

    /* command line parsing */
    char *fn = argv[1];
    long p = atol(argv[2]);
    GEN gen_p = stoi(p);
    long beta = atol(argv[3]);
    beta = 2*(beta/2)+1;

    /* read dimension */
    long ndims, *dims, dim;
    FILE *fp = fopen(fn, "rb");
    fread(&ndims, sizeof(long), 1, fp);
    dims = (long *) malloc(ndims*sizeof(long));
    fread(dims, sizeof(long), ndims, fp);
    dim = dims[0];
    free(dims);

    /* offsets */
    count = 0;
    long *ind = (long *) malloc(dim*dim*sizeof(long));
    for (i = 0; i < dim; i++)
        for (j = i; j < dim; j++)
            ind[dim*i+j] = ind[dim*j+i] = count++;
    long len = count;

    /* read in matrix */
    struct stat filestat;
    stat(fn, &filestat);
    long datacount = (filestat.st_size-(ndims+1)*sizeof(long))/sizeof(long);
    long *data = (long *) malloc(datacount*sizeof(long));
    GEN *T = (GEN *) malloc(len*sizeof(GEN));
    fread(data, sizeof(long), datacount, fp);
    pari_relink(T, data, len);

    /* diagonalize */
    GEN DiagT;
    diagonalize_Zp(T, &DiagT, gen_p, dim, len, ind);

    /* calculate Hnp */
    long pseudodim = 2*((dim+1)/2);

    long *block, maxdeg; GEN *Hnp;
    //calculate_Hnp(gen_p, dim, &Hnp, &block, &maxdeg);
    calculate_Hnp(gen_p, pseudodim, &Hnp, &block, &maxdeg);

    /* initialize */
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* allocate count buffers */
    //long local_subslice_size = 2*(dim+1)*(maxdeg+1);
    long local_subslice_size = 2*(pseudodim+1)*(maxdeg+1);
    long local_slice_size = (beta+1)*local_subslice_size;
    long global_size = nprocs*local_slice_size;

    long *local_term_counts = (long *) calloc(local_slice_size, sizeof(long));
    long *global_term_counts = NULL;
    if (rank == 0)
        global_term_counts = (long *) malloc(global_size*sizeof(long));

    /* iterators */
    long *lambda = (long *) malloc(dim*sizeof(long));
    long *pow = (long *) malloc(dim*sizeof(long));
    long *elem = (long *) malloc(len*sizeof(long));
    long *mod = (long *) malloc(len*sizeof(long));

    /* matrices */
    GEN *G = (GEN *) malloc(len*sizeof(GEN));
    GEN *TInvG = (GEN *) malloc(len*sizeof(GEN));
    GEN *ModPTInvG = (GEN *) malloc(len*sizeof(GEN));
    GEN DiagModPTInvG; long sp, lambdap;
    pari_sp ltop, lbot;
    ltop = avma;

    if (rank == 0) {
        pari_printf("T = \n");
        print_symmetric_matrix(T, dim, ind);
    }

    /* main loop */
    for (i = 0; i <= beta/2; i++) {
        /* split up work */
        long seq_size;
        calculate_work(rank, nprocs, p, dim, i, lambda, elem, &seq_size, ind);

        if (rank == 0)
            pari_printf("seq_size = %ld\n", seq_size);

        /* inner loop */
        count = 0; done = false;
        do {
            for (j = 0; j < dim; j++) {
                pow[j] = 1;
                for (k = 0; k < lambda[j]; k++)
                    pow[j] *= p;
            }
            for (j = 0; j < dim; j++)
                mod[ind[dim*j+j]] = 1;
            for (j = 0; j < dim; j++)
                for (k = j+1; k < dim; k++)
                    mod[ind[dim*j+k]] = pow[k];

            do {
                /* garbage collection */
                if (count >= seq_size) {done = true; break;}
                if (subcount > 1000) {
                    lbot = avma;
                    gerepile(ltop, lbot, NULL);
                    ltop = avma;
                    subcount = 0;
                }

                /* main calculation */
                populate_G(G, dim, pow, elem, ind);
                calculate_TInvG(G, DiagT, dim, TInvG, ind);

                if (is_p_integral(TInvG, gen_p, dim, ind)) {
                    reduce_modp(TInvG, ModPTInvG, gen_p, dim, ind);
                    diagonalize_Fp(ModPTInvG, &DiagModPTInvG, gen_p, dim, len, ind);
                    calculate_p_local_parameters(DiagModPTInvG, gen_p, dim, &sp, &lambdap);
                    for (l = 2*i; l <= beta; l++) {
                        k = l-2*i;
                        if (k <= maxdeg) {
                            //long index = 2*l*(dim+1)*(maxdeg+1)+2*sp*(maxdeg+1)+(maxdeg+1)*lambdap+k;
                            long index = 2*l*(pseudodim+1)*(maxdeg+1)+2*sp*(maxdeg+1)+(maxdeg+1)*lambdap+k;
                            local_term_counts[index]++;
                        }
                    }
                }

                count++; subcount++;
            } while (next_elem(elem, mod, len));
            if (done) break;
            memset(elem, 0, len*sizeof(long));
        } while (next_part(lambda, dim));
    }
        
    MPI_Gather(local_term_counts, local_slice_size, MPI_LONG, 
               global_term_counts, local_slice_size, MPI_LONG,
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        memset(local_term_counts, 0, local_slice_size*sizeof(long));
        long proc;

        //for (i = 0; i <= dim; i++)
        for (i = 0; i <= pseudodim; i++)
            for (j = 0; j <= 1; j++)
                for (k = 0; k <= maxdeg; k++)
                    for (l = 0; l <= beta; l++) {
                        long local_index = local_subslice_size*l+2*i*(maxdeg+1)+(maxdeg+1)*j+k;
                        for (proc = 0; proc < nprocs; proc++) {
                            long global_index = local_slice_size*proc+local_index;
                            local_term_counts[local_index] += global_term_counts[global_index];
                        }
                    }

        GEN S = zerovec(beta+1);

        //for (i = 0; i <= dim; i++)
        /* for (i = 0; i <= pseudodim; i++)
            for (j = 0; j <= 1; j++)
                for (k = 0; k <= maxdeg; k++)
                    for (l = 0; l <= beta; l++) {
                        long local_subslice_index = 2*i*(maxdeg+1)+(maxdeg+1)*j+k;
                        long term_count = local_term_counts[l*local_subslice_size+local_subslice_index];
                        gel(S,l+1) = gadd(gel(S,l+1),gmulgs(Hnp[local_subslice_index],term_count));
                    } */

        pari_printf("\nCalculations\n\n");
        for (l = 0; l <= beta; l++) {
            pari_printf("phi(T;%ld^%ld) = %ld^%ld * (", p, l, p, (l+1)/2);
            for (i = 0; i <= pseudodim; i++)
                for (j = 0; j <= 1; j++)
                    for (k = 0; k <= maxdeg; k++) {
                        long local_subslice_index = 2*i*(maxdeg+1)+(maxdeg+1)*j+k;
                        long term_count = local_term_counts[l*local_subslice_size+local_subslice_index];
                        gel(S,l+1) = gadd(gel(S,l+1),gmulgs(Hnp[local_subslice_index],term_count));
                        if (term_count > 0)
                            pari_printf(" %ld*%Ps +", term_count, Hnp[local_subslice_index]);
                    }
            pari_printf(" )\n");
        }

        for (l = 0; l <= beta; l++)
            gel(S,l+1) = gmul(gel(S,l+1),gpowgs(gen_p,(l+1)/2));

        pari_printf("\nResults\n\n");
        for (l = 0; l <= beta; l++)
            pari_printf("phi(T;%ld^%ld) = %Ps\n", p,l,gel(S,l+1));

        free(global_term_counts);
    }

    /* clean up */
    free(ind); free(data); free(local_term_counts); free(block);
    free(lambda); free(pow); free(elem); free(mod);
    free(T); free(G); free(TInvG);

    /* finalize */
    pari_close(); MPI_Finalize();
    return 0;
}