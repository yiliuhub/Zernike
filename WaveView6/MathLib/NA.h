#ifndef NA_H
#define NA_H

#include <stdio.h>
#include "NAbasicdef.h"


typedef double REAL;

void *vminit         /* create an empty vector/matrix list ...........*/
        (
         void
        );

class NA
{
public:
	NA();
	~NA();
	REAL aPlusb(double a, double b);
	int gauss(int     mod,           /* Modus: 0, 1, 2, 3 ...............*/
           int     n,             /* Dimension of matrix .............*/
           REAL *  mat[],         /* Input matrix ....................*/
           REAL *  lumat[],       /* LU decomposition ................*/
           int     perm[],        /* row remutation vector ...........*/
           REAL    b[],           /* right hand side .................*/
           REAL    x[],           /* solution of the system ..........*/
           int *   signd          /* sign of the permutation .........*/
          );

	int gaudec
           (
            int     n,            /* size of matrix ..................*/
            REAL *  mat[],        /* Input matrix ....................*/
            REAL *  lumat[],      /* matrix decomposition ............*/
            int     perm[],       /* row interchanges ................*/
            int *   signd         /* sign of perm ....................*/
           );

protected:

};

#endif