#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <pari/pari.h>
#include "common.h"




long
count_t_INT (GEN elem)
{ 
    return lg(elem);
}

long
count_t_FRAC (GEN elem)
{
    long count = 3; long j;
    for (j = 1; j < 3; j++)
        count += lg((GEN) elem[j]);
    return count;
}

long
count_t_PADIC (GEN elem)
{
    long count = 5; long j;
    for (j = 2; j < 5; j++)
        count += lg((GEN) elem[j]);
    return count;
}

long 
pari_count(GEN elem)
{
    long count = 0;

    switch (typ(elem)) {
        case t_INT :
            count = count_t_INT(elem);
            break;
        case t_FRAC :
            count = count_t_FRAC(elem);
            break;
        case t_PADIC :
            count = count_t_PADIC(elem);
            break;
    }

    return count;
}

void
pack_t_INT (GEN elem, long **ptr)
{
    long lg = lg(elem);
    memcpy(*ptr, elem, lg*sizeof(long));
    *ptr += lg;
}

void
pack_t_FRAC (GEN elem, long **ptr) 
{
    long *tail = *ptr + 3; long lg, i;
    (*ptr)[0] = elem[0];
    for (i = 1; i < 3; i++) {
        (*ptr)[i] = lg = lg((GEN) elem[i]);
        memcpy(tail, (GEN) elem[i], lg*sizeof(long));
        tail += lg;
    }
    *ptr = tail;
}

void
pack_t_PADIC (GEN elem, long **ptr)
{
    long *tail = *ptr + 5; long lg, i;
    (*ptr)[0] = elem[0]; (*ptr)[1] = elem[1];
    for (i = 2; i < 5; i++) {
        (*ptr)[i] = lg = lg((GEN) elem[i]);
        memcpy(tail, (GEN) elem[i], lg*sizeof(long));
        tail += lg;
    }
    *ptr = tail;
}

void pari_pack (GEN * V, long **data, long *datacount, long len)
{
    long i; *datacount = 0;

    for (i = 0; i < len; i++)
        *datacount += pari_count(V[i]);

    *data = (long *) malloc(*datacount*sizeof(long));
    long *ptr = *data;

    for (i = 0; i < len; i++)
        switch (typ(V[i])) {
            case t_INT :
                pack_t_INT(V[i],&ptr);
                break;
            case t_FRAC :
                pack_t_FRAC(V[i],&ptr);
                break;
            case t_PADIC :
                pack_t_PADIC(V[i],&ptr);
                break;
        }

}

long
relink_t_INT (GEN *elem, long *ptr)
{
    long lg = lg((GEN) ptr);
    *elem = (GEN) ptr;
    return lg;
}

long
relink_t_FRAC (GEN *elem, long *ptr)
{
    long i; long count = 3;
    *elem = (GEN) ptr; long *tail = ptr + 3;
    for (i = 1; i < 3; i++) {
        long lg = ptr[i];
        ptr[i] = (long) tail; 
        tail += lg; count += lg;
    }
    return count;
}

long
relink_t_PADIC (GEN *elem, long *ptr)
{
    long i; long count = 5;
    *elem = (GEN) ptr; long *tail = ptr + 5;
    for (i = 2; i < 5; i++) {
        long lg = ptr[i];
        ptr[i] = (long) tail; 
        tail += lg; count += lg;
    }
    return count;
}

void
pari_relink (GEN *V, long *data, long len)
{
    long i, count;
    for (i = 0; i < len; i++) {
        switch (typ(data)) {
            case t_INT :
                count = relink_t_INT(V,data);
                break;
            case t_FRAC :
                count = relink_t_FRAC(V,data);
                break;
            case t_PADIC :
                count = relink_t_PADIC(V,data);
                break;
        }
        V += 1; data += count;
    }
}

void
diagonalize_Fp(GEN *A, GEN *D, GEN p, long n, long len, long *ind)
{
    pari_sp ltop = avma;
    long i, j, k; long *perm; 

    perm = (long *) pari_malloc(n*sizeof(long));
    for (i = 0; i < n; i++) perm[i] = i;

    GEN *L; L = (GEN *) pari_malloc(len*sizeof(GEN));
    memcpy(L,A,len*sizeof(GEN));
    
    for (i = 0; i < n-1; i++) {
        long rpiv = -1; long cpiv = -1; bool foundnz = false;

        for (j = i; j < n; j++)
            if (!gequal0(L[ind[n*perm[j]+perm[j]]]))
                {rpiv = j; cpiv = j; foundnz = true; break;}

        if (!foundnz)
            for (j = i+1; j < n; j++) {\
                if (foundnz) break; 
                for (k = i; k < j; k++)
                    if (!gequal0(L[ind[n*perm[j]+perm[k]]]))
                        {rpiv = j; cpiv = k; foundnz = true; break;} 
        }

        if (rpiv > -1) {
            if (cpiv != i) {long temp = perm[i]; perm[i] = perm[cpiv]; perm[cpiv] = temp;}
            if (rpiv != cpiv) {
                for (k = i; k < n; k++) 
                    L[ind[n*perm[i]+perm[k]]] = gadd(L[ind[n*perm[i]+perm[k]]],L[ind[n*perm[rpiv]+perm[k]]]);
                L[ind[n*perm[i]+perm[i]]] = gadd(L[ind[n*perm[i]+perm[i]]],L[ind[n*perm[i]+perm[rpiv]]]);
            }

            for (j = i+1; j < n; j++) {
                GEN quot = gdiv(L[ind[n*perm[j]+perm[i]]],L[ind[n*perm[i]+perm[i]]]);
                for (k = i+1; k <= j; k++)
                    L[ind[n*perm[j]+perm[k]]] = gsub(L[ind[n*perm[j]+perm[k]]],gmul(quot,L[ind[n*perm[i]+perm[k]]]));
            }
        }
    }
    
    pari_sp lbot = avma;
    *D = zerovec(n);
    for (i = 0; i < n; i++)
        gel(*D,i+1) = gcopy(L[ind[n*perm[i]+perm[i]]]);
    *D = gerepile(ltop, lbot, *D);
    pari_free(perm); pari_free(L);
}

void
diagonalize_Zp(GEN *A, GEN *D, GEN gen_p, long n, long len, long *ind)
{
    pari_sp ltop = avma;
    long i, j, k; long *perm;

    perm = (long *) pari_malloc(n*sizeof(long));
    for (i = 0; i < n; i++) perm[i] = i;
        
    GEN *L; L = (GEN *) pari_malloc(len*sizeof(GEN));
    memcpy(L,A,len*sizeof(GEN));
    
    for (i = 0; i < n-1; i++) {
        long rpiv = -1; long cpiv = -1; long min;

        for (j = i; j < n; j++)
            if (!gequal0(L[ind[n*perm[j]+perm[j]]]))
                { long val = Q_pval(L[ind[n*perm[j]+perm[j]]],gen_p);
                  if ((rpiv == -1) || (val < min)) 
                  {rpiv = j; cpiv = j; min = val;} }

        for (j = i+1; j < n; j++) {
            for (k = i; k < j; k++)
                if (!gequal0(L[ind[n*perm[j]+perm[k]]]))
                    { long val = Q_pval(L[ind[n*perm[j]+perm[k]]],gen_p);
                      if ((rpiv == -1) || (val < min)) 
                      {rpiv = j; cpiv = k; min = val;} }  
        }

        if (rpiv > -1) {
            if (cpiv != i) {long temp = perm[i]; perm[i] = perm[cpiv]; perm[cpiv] = temp;}
            if (rpiv != cpiv) {
                for (k = i; k < n; k++) 
                    L[ind[n*perm[i]+perm[k]]] = gadd(L[ind[n*perm[i]+perm[k]]],L[ind[n*perm[rpiv]+perm[k]]]);
                L[ind[n*perm[i]+perm[i]]] = gadd(L[ind[n*perm[i]+perm[i]]],L[ind[n*perm[i]+perm[rpiv]]]);
            }

            for (j = i+1; j < n; j++) {
                GEN quot = gdiv(L[ind[n*perm[j]+perm[i]]],L[ind[n*perm[i]+perm[i]]]);
                for (k = i+1; k <= j; k++)
                    L[ind[n*perm[j]+perm[k]]] = gsub(L[ind[n*perm[j]+perm[k]]],gmul(quot,L[ind[n*perm[i]+perm[k]]]));
            }
        }
    }
    
    pari_sp lbot = avma;
    *D = zerovec(n);
    for (i = 0; i < n; i++)
        gel(*D,i+1) = gcopy(L[ind[n*perm[i]+perm[i]]]);
    *D = gerepile(ltop, lbot, *D);
    pari_free(perm); pari_free(L);
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