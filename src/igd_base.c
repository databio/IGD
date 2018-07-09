//=====================================================================================
//Common igd struct, parameters, functions
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"
//-------------------------------------------------------------------------------------

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

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
long double _hypergeo(long long n11, long long n1_, long long n_1, long long n)
{
    //***DEBUG***
    return expl(_lbinom(n1_, n11) + _lbinom(n-n1_, n_1-n11) - _lbinom(n, n_1));
}

// incremental version of hypergenometric distribution
long double _hypergeo_acc(long long n11, long long n1_, long long n_1, long long n, _hgacc_t *aux)
{
    if (n1_ || n_1 || n) {
        aux->n11 = n11; aux->n1_ = n1_; aux->n_1 = n_1; aux->n = n;
    } else { // then only n11 changed; the rest fixed
        if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
            if (n11 == aux->n11 + 1) { // incremental
                aux->p *= (long double)(aux->n1_ - aux->n11) / n11
                    * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
                aux->n11 = n11;
                return aux->p;
            }
            if (n11 == aux->n11 - 1) { // incremental
                aux->p *= (long double)aux->n11 / (aux->n1_ - n11)
                    * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
                aux->n11 = n11;
                return aux->p;
            }
        }
        aux->n11 = n11;
    }
    aux->p = _hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);

    return aux->p;
}

long double _kt_fisher_exact(long long n11,
                             long long n12,
                             long long n21,
                             long long n22,
                             long double *_left,
                             long double *_right,
                             long double *two)
{
    long long i, j, max, min;
    long double p, q, left, right;
    _hgacc_t aux;
    long long n1_, n_1, n;

    n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n

    max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
    min = n1_ + n_1 - n;    // not sure why n11-n22 is used instead of min(n_1,n1_)
    if (min < 0) min = 0; // min n11, for left tail
    *two = *_left = *_right = 1.;

    if (min == max) return 1.; // no need to do test


    q = _hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
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
    if (labs((long) (i - n11)) < labs((long) (j - n11)) && q != 0.0) right = 1. - left + q;
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

