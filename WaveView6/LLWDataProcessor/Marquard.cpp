#include "stdafx.h"

/***************************************************************************
 * File marquardt.c
 *
 * Implemention of Modified Marquardt procedure for function minimization.
 *
 * Adapted from J. C. Nash, "Compact Numerical Methods for Computers"
 * with some modifications
 *
 *     int     choleski_decompos()
 *     void    choleski_back()
 *     double* symmatrix_2_lowtriangle()
 *     double* lowtriangle_2_symmatrix()
 *     double  sumation()
 *     void    derivatives()
 *     int     converge_test()
 *     int     marquardt()
 * 
 * 
 ***************************************************************************/
#include "stdafx.h"

#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "Marquard.h"


/***************************************************************************** 
 *  the max number of parameter so far
 *  circle in 2-D   n = 3 
 *  sphere          n = 4
 *  line   in 3-D   n = 6
 *  circle in 3-D   n = 7
 *  cone            n = 8
 *  cylinder        n = 7
 *
 *  by using MAX_N we can declear the maximum array for vectors and matrix in
 *  marquardt(), so that to elimete malloc() function calls.
 *****************************************************************************/
#define MAX_N  13 
typedef double (*Function)(double*, double*, double*);

/******************************************************************************
 *  choleski_decompos() computes the Choleski decomposition one column at a
 *  time, overwriting the input matrix.  It returns 1 if it succeeds, 0 if it
 *  detects a non-positive diagonal element.  0 indicates that the input
 *  matrix is not (computationally) positive definite.
 *****************************************************************************/
 
int 
choleski_decompos( double *matrix, int row)
{
    int i,k,p;
    int col=row;
 
    for( k=0; k<row; k++ ) 
    {
        for( p=0; p<k; p++ )
            matrix[k*col + k] -= matrix[k*col+p] * matrix[k*col+p];
        if ( matrix[k*col + k] <=0 )
            return( 0 );        /* Not positive definite */
        matrix[k*col + k] = sqrt( matrix[k*col + k] );
        for( i=k+1; i<row; i++ ) {
            for( p=0; p<k; p++ )
                matrix[i*col+k] -= matrix[i*col+p] * matrix[k*col+p];
            matrix[i*col+k]/= matrix[k*col+k];
        }
    }
    return( 1 );
}
 
/****************************************************************************
 *  choleski_back solves a symmetric system given the Choleski factor of the
 *  system.  It assumes the system was (computationally) positive definite.
 *  The solution is carried out via forward elimination followed by back
 *  substitution (modified to use the transpose of the lower triangular
 *  Choleski factor).
 ****************************************************************************/

double* 
choleski_back( double *matrix, double *v, int row)
{
    int i,j;
    int col = row;
 
    /* Forward elimination from lower triangular factor */
 
    for( i=0; i<row; i++) 
    {
        for( j=0; j<i; j++ )
            v[i] -= matrix[i*col+j] * v[j];
        v[i] /= matrix[i*col+i];
    }
 
    /* Back-substitution using transpose of lower triangular factor */
 
    for( i=row-1; i>=0; i--) 
    {
        for( j=i+1; j<row; j++)
            v[i] -= matrix[j*col+i] * v[j];
        v[i] /= matrix[i*col+i];
    }

	return v;
}

/*******************************************************************************
 *   matrix is a n*n array. low_triangle is a n*(n+1)/2 array
 *******************************************************************************/
double*
symmatrix_2_lowtriangle( double *matrix, double *low_triangle, int n)
{
   int i; 
   int j;
 
   for (i = 0; i< n; i++)
   {
       for (j = 0; j <=i; j++)
          low_triangle[i*(i+1)/2+j] = matrix[i*n+j];
   }
   return low_triangle;
}

/******************************************************************************* 
 *   matrix is a n*n array. low_triangle is a n*(n+1)/2 array
 *******************************************************************************/ 
double*
lowtriangle_2_symmatrix( double *low_triangle, double *matrix, int n  )
{
   int i; 
   int j; 
 

   for (i = 0; i< n; i++) 
   { 
       for (j = 0; j <i; j++) 
          matrix[i*n+j] = low_triangle[i*(i+1)/2+j]; 
       for (j = i; j < n; j++)
          matrix[i*n+j] = low_triangle[j*(j+1)/2+i]; 
   } 
   return matrix;
} 


/****************************************************************************
 * void derivatives(points, m, fun_vec, der_fun, b, n, a, v)
 *      Calculates Jacobian Matrix and J'J and vector J'f
 * 
 *   points  ------ are stored as x[i] = points[3*i], y[i] = points[3*i+1]
 *                               z[i] = points[3*i+2]
 *   m       ------ number of points
 *   fun_vec  ------ vector of Fi() for i from 0 to m 
 *   der_fun ------ d(Fi)/d(b[0]) .... d(Fi)/d(b[n])
 *   b       ------ solution vector 
 *   n       ------ number of items in b
 *   a       ------ matrix J'J
 *   v       ------ vector J'f
 ****************************************************************************/
void
derivatives(
   double   *points,
   int      m,
   double   *fun_vec,
   double   (*der_fun)(double*, double*, double*),
   double   *b,
   int      n,
   double   *a,
   double   *v)
{

   int i, j, k, q;
   double s;
   double d[MAX_N];
   double point[3];

   for (j = 0; j < n*(n+1)/2; j++)
       a[j] = 0;
   for (j = 0; j < n; j++)
       v[j] = 0;

   for (i = 0; i< m; i++)
   {
      point[0] = points[3*i];
      point[1] = points[3*i+1];
      point[2] = points[3*i+2];

      der_fun(point, b, d);
/*
      s = fun(point, b);
*/
      s = fun_vec[i];

      for (j=0; j<n; j++)
      {
         v[j]+=d[j]*s; /* accumulate J'f */
         q = j*(j+1)/2;
         for (k= 0; k<=j; k++)
             a[q+k]+=d[j]*d[k]; /* accumulate J'J */

      }
   }
}

/****************************************************************************
 * double sumation(points, m, fun, b, fun_vec)
 *   calculates the sum of the distance squares from points to the object
 *
 *   points  ------ are stored as x[i] = points[3*i], y[i] = points[3*i+1]
 *                               z[i] = points[3*i+2]
 *   fun     ------ Fi()*Fi() at point i
 *   m       ------ number of points
 *   b       ------ solution vector so far
 *   fun_vec   ------ vector to save Fi
 ****************************************************************************/

double
sumation(
   double *points,
   int m,
   double (*fun)(double*, double*),
   double* b,
   double *fun_vec)
{
 
   int i;
   double fi;
   double sum=0.0;
   double point[3];
 
   for (i = 0; i< m; i++)
   {
      point[0] = points[3*i];
      point[1] = points[3*i+1];
      point[2] = points[3*i+2];
      fi= fun(point,b);
      fun_vec[i] = fi;
      sum +=fi*fi;

   }   
   return sum;
}

/**************************************************************************
 *  Test if (new_v - old_v) has max norm <= tolerance (tolerance >= 0) 
 *  or if   component-wise, (new_v-old_v) <= -tolerance*old_v (tolerance < 0).
 **************************************************************************/

int
converge_test(
   double* new_v,
   double* old_v,
   int n,
   double tolerance)
{
   int i;

   if (tolerance <0.0)
   /* Relative change test; use -tolerance as convergence factor */
   {
      for (i = 0; i < n; i++)
          if (fabs(new_v[i]-old_v[i]) > -tolerance*fabs(old_v[i]))
             return 0;
   }
   else
   {
      /* Absolute change test */
      for (i = 0; i < n; i++)
          if (fabs(new_v[i]-old_v[i]) > tolerance)
             return 0;
   }
   return 1;
}

/****************************************************************************
 *  Modified Marquardt procedure for function minimization.
 *
 *  Adapted from J. C. Nash, Compact Numerical Methods for Computers:
 *  Linear Algebra and Function Minimisation, Adam Hilger Ltd., Bristol,
 *  England, 1979.
 *
 *  Explanation of parameters to marquardt():
 *      points  ------ points used in evaluation of sum and derivs
 *      fun     ------ Fi() at point i
 *      der_fun ------ d(Fi)/d(b[0]) .... d(Fi)/d(b[n])
 *      m       ------ number of points
 *      b       ------ input: the starting guess; output the solution vector
 *      n       ------ number of items in b
 *      tolerance   ------ user-supplied convergence factor ( <0 => relative)
 *      iterateLimit  ------ user-supplied limit on function (sum) evaluations
 *      sum_sq  ------ summation of distance square
 *
 *  The return value for the function is a ResultCode indicating success or
 *  the reason for failure:
 *      1   - success
 *      2  - limit reached on function evaluation limit
 *      0 - unable to allocate working memory from the heap
 ****************************************************************************/
int
marquardt(
   double   *points,
   int      m, /*# of residuals--- number of points */
   double (*fun)(double*, double*),
   double ( *der_fun)(double*, double*, double*),
   double*  b, /* initial value */
   int      n, /* number of parameters */
   double   tolerance,
   int      iterateLimit,
   double   *sum_sq)
{
   double inc=7.0;   // = 10.0;
   double dec=0.618; // = 0.4;
   double lambda=0.0001;
   double shrink_ratio = 3;
   double lambdaLimit = 1.e-030;
   double phi=1;
   int ifn = 0;
   int ig = 0;

   double p;
   double p0 = 0.0;
   double a[MAX_N*(MAX_N+1)/2];   /* n*(n+1)/2 is needed */
   double a2[MAX_N*MAX_N];       /* n*n is needed */
   double c[MAX_N*(MAX_N+1)/2];  /* n*(n+1)/2 is needed */
   double v[MAX_N];              /* n is needed */
   double x[MAX_N];              /* n is needed */
   double d[MAX_N];              /* n is needed */
   double *fun_vec;

 

   int converged = 0;
   int progressed ;
   int i, j, q;
   int count;

   if(tolerance > 0.00000000001)
	   tolerance = 0.00000000001;
   else if(tolerance <= 1.e-15)
	   tolerance = 0.00000000001;


   fun_vec = (double*) malloc(sizeof(double)*m);
   if (fun_vec == NULL)
      return 0;

   p = sumation(points, m, fun, b, fun_vec);
   ifn++;


   /*2*/
   while (!converged &&(iterateLimit==0||ifn <iterateLimit))
   {
       ig++;
	   p0=p;

	   if(fabs(lambda) > lambdaLimit)
		   lambda*=dec;
		

       /* now compute J'J, J'f */

       /*3*/  
       for (j = 0; j < n*(n+1)/2; j++)
          a[j] = 0;
       for (j = 0; j < n; j++)
          v[j] = 0;

       /* 4 */
       derivatives(points, m, fun_vec, der_fun, b, n, a, v);

       for (j = 0; j < n*(n+1)/2; j++)
           c[j] = a[j];

       for (j = 0; j < n; j++)
           d[j] = b[j];
       progressed = 0;
       while (!progressed &&(iterateLimit==0||ifn <iterateLimit))  
       {
		   ifn++;
           for (j=0; j<n; j++)
           {
               q=j*(j+1)/2+j;
               a[q]=c[q]*(1+lambda)+phi*lambda;
               x[j]=-v[j];
           }

           lowtriangle_2_symmatrix(a, a2, n);
           if (choleski_decompos(a2, n)) /* positive defined */
           {
              choleski_back(a2, x, n);
              count=0;

              for (i = 0; i < n; i++)
              {
                  b[i]=d[i]+x[i];
                  if (b[i]==d[i])
                     count = count+1;
              }
              p = sumation(points, m, fun, b, fun_vec);
              
              progressed=p<p0;
              if ( progressed || p == p0 ) 
              {
                  converged = converge_test( b, d,n, tolerance );
                  if (converged)
                  {
					  if(p0 - p < 1.0e-15)
						  break;
					  else
						  converged = 0;

					  p0 = p;
                  }
              }
           }

           if (!progressed)
		   {
				lambda*=inc;
				if(lambda > 1.0e6 && p == p0)
				{
					lambda = 0.0001;
					inc = 1 + (inc - 1)/shrink_ratio;
				}
		   }

       } 

       if (!progressed) 
          for (i = 0; i < n; i++)
            b[i]=d[i];

   }
   if (sum_sq !=NULL)
      *sum_sq = p0;

   free(fun_vec);

   if (converged)
      return 1;
   else
      return 2; /* limit reached */
}

