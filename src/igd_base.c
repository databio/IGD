//=====================================================================================
//Common igd struct, parameters, functions
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"

int compare_iidx(const void *a, const void *b)
{
    struct igd_data *pa = (struct igd_data *) a;
    struct igd_data *pb = (struct igd_data *) b;
    return pa->i_idx - pb->i_idx;
}

int compare_rend(const void *a, const void *b)
{
    struct igd_data *pa = (struct igd_data *) a;
    struct igd_data *pb = (struct igd_data *) b;
    return pa->r_end - pb->r_end;
}

int compare_iidx1(const void *a, const void *b)
{
    struct igd_data1 *pa = (struct igd_data1 *) a;
    struct igd_data1 *pb = (struct igd_data1 *) b;
    return pa->i_idx - pb->i_idx;
}

int compare_rend1(const void *a, const void *b)
{
    struct igd_data1 *pa = (struct igd_data1 *) a;
    struct igd_data1 *pb = (struct igd_data1 *) b;
    return pa->r_end - pb->r_end;
}

int compare_qidx(const void *a, const void *b)
{
    struct query_data *pa = (struct query_data *) a;
    struct query_data *pb = (struct query_data *) b;
    return pa->q_idx - pb->q_idx;
}

int compare_rstart(const void *a, const void *b)
{
    struct query_data *pa = (struct query_data *) a;
    struct query_data *pb = (struct query_data *) b;
    return pa->r_start - pb->r_start;
}

int compare_rstart2(const void *a, const void *b)
{
    struct igd_data2 *pa = (struct igd_data2 *) a;
    struct igd_data2 *pb = (struct igd_data2 *) b;
    return pa->r_start - pb->r_start;
}

int compare_midx(const void *a, const void *b)
{
    struct igd_mix *pa = (struct igd_mix *) a;
    struct igd_mix *pb = (struct igd_mix *) b;
    return pa->m_idx - pb->m_idx;
}

char** str_split_t( char* str, int nItems)
{
    char **splits;
    char *tmp;
    int i;
    if (str == NULL)
        return NULL;
    else {    
        splits = malloc((nItems+1)*sizeof(*splits)); 
        i=0;
        tmp = strtok(str, "\t");
        while(tmp!=NULL && i<nItems){
            splits[i] = tmp;
            tmp = strtok(NULL, "\t");
            i++;
        }
    }
    //printf("%s %s %s \n", splits[0], splits[1], splits[2]);
    return splits;
}

void str_splits( char* str, int *nmax, char **splits)
{   //tsv 
    splits[*nmax] = NULL;
    splits[0] = str;
    char *ch = str;    
    int ns = 1;    
    do {
        if (*ch == '\t'){
            splits[ns++] = &ch[1];
            *ch = '\0';
        }
        ch++;
    } while (*ch != '\0' && ns < *nmax+1);
    *nmax = ns;
}

char** str_split( char* str, char delim, int *nmax)
{   //slightly faster than _t
    char** splits;
    char* ch;    
    int ns;
    if (str == NULL || delim=='\0')
        return NULL;
    else {
        splits = malloc((*nmax+1) * sizeof(*splits));
        splits[*nmax] = NULL;
        ch = str;
        ns = 1;
        splits[0] = str;
        do {
            if (*ch == delim)
            {
                splits[ns++] = &ch[1];
                *ch = '\0';
            }
            ch++;
        } while (*ch != '\0' && ns < *nmax+1);
    }
    *nmax = ns;
    return splits;
}

//-------------------------------------------------------------------------------------
//This section is taken from giggle package
/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 * kfunc.h
 *  Created on: May 1, 2015
 *      Author: nek3d
 */

// log\binom{n}{k}
long double _lbinom(long long n, long long k)
{
    if (k == 0 || n == k) return 0;
    return lgammal(n+1) - lgammal(k+1) - lgammal(n-k+1);
}

// nA   nC   | nAC
// nB   nD   | nBD
//-----------+----
// nAB  nCD  | n

// hypergeometric distribution
long double _hypergeo(long long nA, long long nAC, long long nAB, long long n)
{   //n:   population size 
    //nAB: number of draws
    //nA:  observed success 
    //nAC: success states in n
    return expl(_lbinom(nAC, nA) + _lbinom(n-nAC, nAB-nA) - _lbinom(n, nAB));
}

// incremental version of hypergenometric distribution
long double _hypergeo_acc(long long nA, long long nAC, long long nAB, long long n, _hgacc_t *aux)
{
    if (nAC || nAB || n) {
        aux->nA = nA; aux->nAC = nAC; aux->nAB = nAB; aux->n = n;
    } else { // then only nA changed; the rest fixed
        if (nA%11 && nA + aux->n - aux->nAC - aux->nAB) {
            if (nA == aux->nA + 1) { // incremental
                aux->p *= (long double)(aux->nAC - aux->nA) / nA
                    * (aux->nAB - aux->nA) / (nA + aux->n - aux->nAC - aux->nAB);
                aux->nA = nA;
                return aux->p;
            }
            if (nA == aux->nA - 1) { // incremental
                aux->p *= (long double)aux->nA / (aux->nAC - nA)
                    * (aux->nA + aux->n - aux->nAC - aux->nAB) / (aux->nAB - nA);
                aux->nA = nA;
                return aux->p;
            }
        }
        aux->nA = nA;
    }
    aux->p = _hypergeo(aux->nA, aux->nAC, aux->nAB, aux->n);

    return aux->p;
}

long double _kt_fisher_exact(long long nA,
                             long long nC,
                             long long nB,
                             long long nD,
                             long double *_left,
                             long double *_right,
                             long double *two)
{
    long long i, j, max, min;
    long double p, q, left, right;
    _hgacc_t aux;
    long long nAC, nAB, n;

    nAC = nA + nC; nAB = nA + nB; n = nA + nC + nB + nD; // calculate nAC, nAB and n

    max = (nAB < nAC) ? nAB : nAC; // max nA, for right tail
    min = nAC + nAB - n;    // not sure why nA-nD is used instead of min(nAB,nAC)
    if (min < 0) min = 0; // min nA, for left tail
    *two = *_left = *_right = 1.;

    if (min == max) return 1.; // no need to do test


    q = _hypergeo_acc(nA, nAC, nAB, n, &aux); // the probability of the current table
    if (q < 1e-200) q = 1e-200;

    // left tail
    p = _hypergeo_acc(min, 0, 0, 0, &aux);
    for (left = 0., i = min + 1; p < 0.99999999 * q && i<=max; ++i) // loop until underflow
        left += p, p = _hypergeo_acc(i, 0, 0, 0, &aux);
    --i;
    if (p < 1.00000001 * q) left += p;
    else --i;
    // right tail
    p = _hypergeo_acc(max, 0, 0, 0, &aux);
    for (right = 0., j = max - 1; p < 0.99999999 * q && j>=0; --j) // loop until underflow
        right += p, p = _hypergeo_acc(j, 0, 0, 0, &aux);
    ++j;
    if (p < 1.00000001 * q) right += p;
    else ++j;
    // two-tail
    *two = left + right;
    if (*two > 1.) *two = 1.;
    // adjust left and right
    if (labs((long) (i - nA)) < labs((long) (j - nA)) && q != 0.0) right = 1. - left + q;
    else left = 1.0 - right + q;
    *_left = left; *_right = right;
    return q;
}

double log2fc(double ratio)
{
    if (fabs(ratio) < 0.0001)
        return 0.0;

    if (ratio < 1) {
        ratio = 1.0/ratio;
        return -1.0 * log2(ratio);
    }

    return log2(ratio);
}

long double neglog10p(long double sig)
{
    if (fabsl(sig) < -DBL_MAX)
        return 10.0;
    return -1.0 * log10l(sig);
}

