/*
 TO COMPILE USE THE CODE:
 R CMD SHLIB Cal_nji.c -lgsl -lgslcblas
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_eigen.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "R.h"
#include "Rmath.h"

/* */
void Cal_nji_loop(int *nn,
                  int *nc,
                  int *Jval,
                  double typeVec[],
                  double marksVec[],
                  double xVec[],
                  double yVec[],
                  double ulxExt[],
                  double llxExt[],
                  double ulyExt[],
                  double llyExt[],
                  double nji[])
{
    int i, j, k;
    
    gsl_matrix *n_ji = gsl_matrix_calloc(*nc, *Jval);
    
    for(i = 0; i < *nn; i++)
    {
        j = marksVec[i]-1;
        for(k = 0; k < *nc; k++)
        {
            if(xVec[i]<=ulxExt[k] && xVec[i]>llxExt[k] && yVec[i]<=ulyExt[k] && yVec[i]>llyExt[k])
            {
                gsl_matrix_set(n_ji, k, j, gsl_matrix_get(n_ji, k, j)+1);
            }
        }
        if( ( (i+1) % 1000000 ) == 0)
        {
            Rprintf("Allocating %d out of %d observations \n", i+1, *nn);
            R_FlushConsole();
            R_ProcessEvents();
        }
    }
    
    for(i = 0; i < *Jval; i++)
    {
        for(j = 0; j < *nc; j++)
        {
            nji[i * (*nc) + j] = (int) gsl_matrix_get(n_ji, j, i);
        }
    }
    
    return;
}







