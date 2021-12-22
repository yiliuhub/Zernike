#include "stdafx.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>


#include "examples/img/CImg_demo.h"
#include "CImageTest.h"
//#include "..\LLWDataProcessor\Marquard.h"

//#include "numerical.h"


using namespace cimg_library;
using namespace std;
bool checkFile(string filename)
{
    //open file, if it opens successfull it exists, if it fails it doesnt exist
    ifstream testFile;
    testFile.open(filename.c_str(),ios_base::binary );
 
    bool bExists = false;
    if(testFile)
    {
        bExists = true;
    }
    testFile.close();
 
    return bExists;
}

int **readFile(ifstream inputFile, int& size)
{
    //local variables
    int temp;
//    size=256;
 
    //dynamic 2d array of size rows x size columns
    int ** p = new int*[size];
   
    for(int i=0; i < size; i++)
        p[i] = new int[size];
 
 
    //read from file into temp student
    for(int j=0;j<size; j++)
    {
        for (int i=0;inputFile>>temp && i<size; i++)
            p[i][j]=temp;
    }   
 
    return p;
}


bool openFile(std::string& fileName)
{
  bool fileOpened = checkFile( fileName );

  return fileOpened;
}
struct metaballs3d {
  float cx1, cy1, cz1, cx2, cy2, cz2, cx3, cy3, cz3;
  inline float operator()(const float x, const float y, const float z) const {
    const float
      x1 = x - cx1, y1 = y - cy1, z1 = z - cz1,
      x2 = x - cx2, y2 = y - cy2, z2 = z - cz2,
      x3 = x - cx3, y3 = y - cy3, z3 = z - cz3,
      r1 = 0.3f*(x1*x1 + y1*y1 + z1*z1),
      r2 = 0.4f*(x2*x2 + y2*y2 + z2*z2),
      r3 = 0.5f*(x3*x3 + y3*y3 + z3*z3);
    float potential = 0;
    if (r1<1.3f) potential+= 1.0f - r1*(r1*(4*r1+17)-22)/9;
    if (r2<1.3f) potential+= 1.0f - r2*(r2*(4*r2+17)-22)/9;
    if (r3<1.3f) potential+= 1.0f - r3*(r3*(4*r3+17)-22)/9;
    return potential;
  }
};


void ConvertImageToSpatialPoint(double* imgBuff, int colSize, int rowSize, Point* pt)
{
    long count = 0;
    for(long i = 0; i < colSize; i++)
    {
        for(long j = 0; j < rowSize; j++)
        {
            pt[count].x = i;
            pt[count].y = j;
            pt[count].z = imgBuff[count];

            count++;
        }
    }
}

void Test_GaussianBestfit(double* imgBuff, int colSize, int rowSize)
{
    int pointSize = colSize*rowSize;
    int useInit = 0;
    int iterateLimit = 1000;
    double tolerance = 0.000001;
    double *extraInfo = NULL;

    Point *points = new Point[colSize*rowSize];
    GaussianAffine2D *para = new GaussianAffine2D[1];

    para->amp    = 250;
    para->ctrX    = 6.5;
    para->ctrY    = 6.5;
    para->a        = -0.2;
    para->b        =  0.0;
    para->c        = -0.2;


    int result = 0;

    if(points != NULL)
    {
        ConvertImageToSpatialPoint(imgBuff, colSize, rowSize, points);

        SimulateImageData(colSize, rowSize, points,
                       para->amp, para->ctrX, para->ctrY, para->a, para->b, para->c);

        result = gaussianAffine2D_bestfit(points, pointSize, tolerance, iterateLimit,
                                                para, useInit, extraInfo);

    }

    delete[] para;
}

void openFile(ifstream& inputFile)
{
    //local variable
    string filename;
 
    //user input
    cout<<"Please enter the file path for the .raw image : ";
    cin>>filename;
    cout<<endl<<endl;
 
    //continue to ask for file if file does not exist
    while (checkFile(filename))
    {
        cout<<"This File does not exist, please enter another filename: ";
        cin>>filename;
        cout<<endl<<endl;
    }
 
    inputFile.open(filename.c_str(),ios_base::binary);
}





//bool readRawImageFile(ifstream &inputFile, vector<vector<unsigned char>> &rawdata)
//{
//   unsigned char temp;
//
// // get length of file:
//  inputFile.seekg (0, ios::end);
//  int length = inputFile.tellg();
//  inputFile.seekg (0, ios::beg);
//
//  // determine the number of rows and cols
//  int size = (int)sqrt((double)length);
//  // check if it is a perfect square
//  if (size*size != length)
//       return false;
//
//  // resize the vector
//  rawdata.resize(size);
//  for (int i=0; i<size; i++)
//       rawdata[i].resize(size);
//
//   //read from file into vector
//   for(int j=0;j<size; j++)
//   {
//      for (int i=0; i<size; i++)
//      {
//          inputFile.get(temp);
//          p[j][i] = temp;
//      }
//   }   
//   return true;
//}

template<class ImageType>
bool readRawImage(const char* fileName, int row, int col, bool bColor, ImageType* buff)
{
    FILE *input;
    ImageType get_char;
    input = fopen(fileName, "r");
    if(!input)
    {
        return false;
    }

    long size = row*col;
    if(bColor)
    {
        size *= 3;
    }

    long pos = 0;

    while((get_char=fgetc(input))!= EOF && pos < size)
    {
        buff[pos++] = (ImageType) get_char;
    }
    fclose(input);

    if( pos != size )
        return false;
    else
        return true;
}


template<class ImageType>
bool readBWImageFromInterlacedRawImage(const char* fileName, int row, int col, bool bColor, ImageType* buff)
{
    FILE *input;
    ImageType get_char;
    input = fopen(fileName, "rb");
   
    ifstream infile;
    infile.open(fileName);

    infile.seekg (0, ios::end);
    int length = infile.tellg();
    infile.seekg (0, ios::beg);


    if(!input)
    {
        return false;
    }

    long size = row*col;
   

    long pos = 0;
    long pos1 = 0;


    char temp;
    for( long i = 0; i < size; i++)
    {
        infile.get( temp );
        infile.get( temp );
        infile.get( temp );

        buff[i] = temp;

        if( temp < 0 )
        {
            buff[i] += 256;
        }

        pos++;
    }

    //pos = 0;
    //size = 3*row*col;

    //while((get_char=infile.get())!= EOF && pos < size)
    //{
    //    if(pos%3 == 0)
    //    {
    //        buff[(long)(pos/3)] = (ImageType) get_char;
    //    }
    //    pos++;
    //}

    fclose(input);

    if( pos != size )
        return false;
    else
        return true;
}





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


double gaussianAffine2D(double amf, double ctrX, double ctrY, double a, double b, double c, double x, double y)
{
    return amf * exp ( a*(x-ctrX)*(x-ctrX) + b*(x-ctrX)*(y-ctrY) + c*(y-ctrY)*(y-ctrY) );
}


double gaussianAffine2D_sum_fun(double* point, double* b)
{
    double x = point[0], y = point[1], z = point[2];

    return z - gaussianAffine2D(b[0], b[1], b[2], b[3], b[4], b[5], x, y);
}
double point_to_gaussianAffine2D_sq_D1(double *param, double x)
{
    double dist = 0.0;
    return dist;
}

double point_to_gaussianAffine2D_sq_D2(double *param, double x)
{
    double dist = 0.0;
    return dist;
}

int initial_gaussianAffine2D(Point* point, int n, GaussianAffine2D* init_gaussian)
{
    long currPixel = 0;
    double    maxX = point[0].x,
            maxY = point[0].y,
            maxZ = point[0].z;

    for(long i = 1; i < n; i++)
    {
        if( point[i].z > maxZ )
        {
            maxX = point[i].x;
            maxY = point[i].y;
            maxZ = point[i].z;
        }
    }

    if( init_gaussian != NULL)
    {
        init_gaussian->amp    = maxZ;
        init_gaussian->ctrX = maxX;
        init_gaussian->ctrY = maxY;
        init_gaussian->a = -1.0;
        init_gaussian->b =  0.0;
        init_gaussian->c = -1.0;
    }

    double sum = 0.0;
    long count = 0;
    for(long i = 0; i < n; i++)
    {
        if(point[i].z < maxZ)
        {
            sum = log(fabs(point[i].z) / maxZ) /
                ((point[i].x - maxX)*(point[i].x - maxX) + (point[i].y - maxY)*(point[i].y - maxY));

            count++;
        }
    }

    init_gaussian->a = init_gaussian->c = sum / (double)count;

    return 1;
}


int gaussianAffine2D_bestfit(Point *points, int n, double tolerance, int iterateLimit,
                        GaussianAffine2D *init_gaussianAffine2D, int useInit, double *extraInfo)
{
    double solVec[6];
   int success;
 
   if (n < 4)
   {
      return FAILURE;
   }


   if(!useInit)
   {
        success = initial_gaussianAffine2D(points, n, init_gaussianAffine2D);

        if (!success)
            return FAILURE;
   }

  
   solVec[0] = init_gaussianAffine2D->amp;
   solVec[1] = init_gaussianAffine2D->ctrX;
   solVec[2] = init_gaussianAffine2D->ctrY;
   solVec[3] = init_gaussianAffine2D->a;
   solVec[4] = init_gaussianAffine2D->b;
   solVec[5] = init_gaussianAffine2D->c;

   success = marquardt((double*)points, n, gaussianAffine2D_sum_fun, gaussianAffine2D_deriv_fun,
                       solVec, 6, tolerance, iterateLimit, extraInfo);

   if (!success)
      return FAILURE;

   init_gaussianAffine2D->amp  = solVec[0];
   init_gaussianAffine2D->ctrX = solVec[1];
   init_gaussianAffine2D->ctrY = solVec[2];
   init_gaussianAffine2D->a = solVec[3];
   init_gaussianAffine2D->b = solVec[4];
   init_gaussianAffine2D->c = solVec[5];
  
   return success;
}


double gaussianAffine2D_dist_fun(double* point, double* para)
{
    double  A = para[0],
            x0 = para[1],
            y0 = para[2],
            a = para[3],
            b = para[4],
            c = para[5];
    double x = point[0], y = point[1], z = point[2];

    double G = gaussianAffine2D(A, x0, y0, a, b, c, x, y);
    return fabs(z - G);
}

/*----------------------------------------------------------------------------------*
 | 2D GAUSSIAN FUNCTION                                                                |
 | Type 2:
 | G(x, y) = A * exp( a*(x-x0)^2 + b*(x-x0)*(y-y0) + c*(y-y0)^2 )                    |
 |                                                                                    |
 *----------------------------------------------------------------------------------*/
double gaussianAffine2D_deriv_fun(double* point, double* para, double* d)
{
    double  A = para[0],
            x0 = para[1],
            y0 = para[2],
            a = para[3],
            b = para[4],
            c = para[5];
    double x = point[0], y = point[1], z = point[2];

    double G = gaussianAffine2D(A, x0, y0, a, b, c, x, y);

    double sign = 0.0;        // sign of G - z

    sign = ( G < z)? 1.0 : -1.0;

    d[0] = sign*(A > 0.0)? G/A : 0.0;                // dG/dA
    d[1] = sign*(2*a*(x-x0)+b*(y-y0))*G;            // dG/dx0
    d[2] = sign*(2*c*(y-y0)+b*(x-x0))*G;            // dG/dy0
    d[3] = sign*(x-x0)*(x-x0)*G;                    // dG/da
    d[4] = sign*(x-x0)*(y-y0)*G;                    // dG/db
    d[5] = sign*(y-y0)*(y-y0)*G;                    // dG/dc

    return 1.0;
}


int gaussian_test_fun()
{
    int result = 1;
    return result;
}



void SimulateImageData(int colSize, int rowSize, Point* imgPt,
                       double amp, double ctrX, double ctrY, double a, double b, double c)
{
    long count = 0;
    for( long row = 0; row < rowSize; row++)
    {
        for(long col = 0; col < colSize; col++)
        {
            imgPt[count].x = col;
            imgPt[count].y = row;
            imgPt[count].z = 0.1*cos( (fabs(col-ctrX)+ fabs(row-ctrY))*PIE/colSize) +
                            gaussianAffine2D(amp, ctrX, ctrY, a, b, c, imgPt[count].x, imgPt[count].y );

            count++;
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//  return 0 -- local maximum intensity is less than threshold
//         1 -- there is only one maximum pixel
//         2 -- there are multiple maximum pixels
//
/////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
int SearchLocalExtrem(CImg<T>& img, int centerX, int centerY, int radius,
                       int& loclMaxX, int& loclMaxY, float threshold = 0, bool bMaximum=true)
{
  int result = 1;

  int w = img.width();
  int h = img.height();

  int x[2] = {centerX - radius, centerX + radius};
  int y[2] = {centerY - radius, centerY + radius};

  if( radius < 0 || x[0]<0 || x[1]>w-1 || y[0]<0 || y[1]>h-1  )
  {
    return false;
  }

  T globalMin, globalMax;
  globalMin = img.min_max( globalMax );

  long maxIndex;
  globalMin = img.min_max_index(globalMax, maxIndex);

  return Get_Local_Max(img, x[0], y[0], 0, x[1], y[1], loclMaxX, loclMaxY, threshold);
}

//
// this function works right only when there is no missing spot(s) near image center
//
template<typename T>
bool FindSteSize ( CImg<T>& img, int& maxPixelX, int& maxPixelY, int& stepX, int& stepY)
{
  bool bOk = true;
  int deltaX = 10;
  int deltaY = 10;

  int x[2], y[2];
  x[0] = (int)( 3*img.width()/8.0 );
  y[0] = (int)( 3*img.height()/8.0 );
  x[1] = (int)( 5*img.width()/8.0 );
  y[1] = (int)( 5*img.height()/8.0 );

  float globalMin, globalMax;
  globalMin = img.min_max( globalMax );

  int count = 0;
  int localMaxX, localMaxY;
   
  int found = Get_Local_Max(img, x[0], y[0], 1, x[1], y[1], maxPixelX, maxPixelY, 0);

  if(found < 0)  {
    return false;
  }

  // stepX
  x[0] = maxPixelX;
  x[1] = x[0] + 4*deltaX;
  y[0] = maxPixelY - deltaY;
  y[1] = maxPixelY + deltaY;

  bOk = false;
  found = 0;
  while( found < 1 && x[0] < img.width() - 4*deltaX)
  {
    x[0] += deltaX;
    x[1] = x[0] + 4*deltaX;

    found = Get_Local_Max(img, x[0], y[0], 1, x[1], y[1], localMaxX, localMaxY, 0.5*globalMax);

    if( (localMaxX - x[0] < deltaX) ||
        (x[1] - localMaxX < deltaX) )
    {
      found = 0;
    }
  }

  if(found > 0)
  {
    stepX = localMaxX - maxPixelX;
    bOk = true;
  }
  else
    return bOk;

  // stepY
  x[0] = maxPixelX - deltaX;
  x[1] = x[0] + deltaX;

  y[0] = maxPixelY;
  y[1] = y[0] + 4*deltaY;

  found = 0;
  bOk = true;
  while(found <1 && y[0] < img.height() - 4*deltaY)
  {
    y[0] += deltaX;
    y[1] = y[0] + 4*deltaX;

    found = Get_Local_Max(img, x[0], y[0], 1, x[1], y[1], localMaxX, localMaxY, 0.5*globalMax);

    if( (localMaxY - y[0] < deltaY) || (y[1] - localMaxY < deltaY) )
    {
      bOk = false;
    }
  }

  if(found > 0)
  {
    stepY = localMaxY - maxPixelY;
    bOk = true;
  }

  return bOk;
}

template<typename T>
void DirectionalLocalExtremSearchX( CImg<T>& img, CImg<T>& resultMatrix, const float threshold,
                                    const int startX, const int startY,
                                    const int gridIndexX, const int gridIndexY,
                                    const unsigned int gridStepX, const unsigned int gridStepY)
{
  int Radius = 2;
  unsigned int Minimum_Grid = 10;

  if( (gridStepX < Minimum_Grid) || (gridStepY < Minimum_Grid) )
  {
      return;
  }


  for(int i = -1; i<2; i+=2)// i = -1 first, then +1
  {
    int centerX = startX;
    int centerY = startY;
    int localMaxX, localMaxY;
    unsigned int indX = gridIndexX, indY = gridIndexY;

    while(centerX>Radius && centerX < img.width()-Radius &&
          centerY>Radius && centerY < img.height()-Radius)
    {
        int count = SearchLocalExtrem(img, centerX, centerY, Radius, localMaxX, localMaxY, threshold);

        if(count > 0)
        {
          int x = localMaxX;
          int y = localMaxY;
          float sum3x3 = img(x-1, y-1) + img(x, y-1) + img(x+1, y-1) +
                         img(x-1, y  ) + img(x, y  ) + img(x+1, y  ) +
                         img(x-1, y+1) + img(x, y+1) + img(x+1, y+1);
          float avgX =  (x*sum3x3  - (img(x-1, y-1) + img(x-1, y)+img(x-1, y+1))
                                   + (img(x+1, y-1) + img(x+1, y)+img(x+1, y+1))
                        )/sum3x3;
           float avgY = (y*sum3x3  - (img(x-1, y-1) + img(x, y-1)+img(x+1, y-1))
                                   + (img(x-1, y+1) + img(x, y+1)+img(x+1, y+1))
                        )/sum3x3;
			float modulevalue = 100*((indX+indY)%2);
			resultMatrix(indX, indY, 0, 0) = modulevalue;
            resultMatrix(indX, indY, 0, 0) += 2000;
            resultMatrix(indX, indY, 0, 1) = avgX;
            resultMatrix(indX, indY, 0, 2) = avgY;
           
            centerX = localMaxX + i*gridStepX;
        }
        else
        {
          centerX += i*gridStepX;
        }

        indX += i;

        if(indX < 0 || indX > resultMatrix.width()-1)
        {
          break;
        }
    }
  }// end of for loop
}

template<typename T>
void DirectionalLocalExtremSearchY( CImg<T>& img, CImg<T>& resultMatrix, const float threshold,
                                    const int startX, const int startY,
                                    const int gridIndexX, const int gridIndexY,
                                    const unsigned int gridStepX, const unsigned int gridStepY)
{
  int Radius = 2;
  unsigned int Minimum_Grid = 10;

  if( (gridStepX < Minimum_Grid) || (gridStepY < Minimum_Grid) )
  {
      return;
  }


  for(int i = -1; i<2; i+=2)// i = -1 first, then +1
  {
    int centerX = startX;
    int centerY = startY;
    int localMaxX, localMaxY;
    unsigned int indX = gridIndexX, indY = gridIndexY;

    while(centerX>Radius && centerX < img.width()-Radius &&
          centerY>Radius && centerY < img.height()-Radius)
    {
        int count = SearchLocalExtrem(img, centerX, centerY, Radius, localMaxX, localMaxY, threshold);

        if(count > 0)
        {
          int x = localMaxX;
          int y = localMaxY;
          float sum3x3 = img(x-1, y-1) + img(x, y-1) + img(x+1, y-1) +
                         img(x-1, y  ) + img(x, y  ) + img(x+1, y  ) +
                         img(x-1, y+1) + img(x, y+1) + img(x+1, y+1);
          float avgX =  (x*sum3x3  - (img(x-1, y-1) + img(x-1, y)+img(x-1, y+1))
                                   + (img(x+1, y-1) + img(x+1, y)+img(x+1, y+1))
                        )/sum3x3;
           float avgY = (y*sum3x3  - (img(x-1, y-1) + img(x, y-1)+img(x+1, y-1))
                                   + (img(x-1, y+1) + img(x, y+1)+img(x+1, y+1))
                        )/sum3x3;
            resultMatrix(indX, indY, 0, 0) = 3000;
            resultMatrix(indX, indY, 0, 1) = avgX;
            resultMatrix(indX, indY, 0, 2) = avgY;
           
            centerY = localMaxY + i*gridStepY;
        }
        else
        {
          centerY += i*gridStepY;
        }

        indY += i;

        if(indY < 0 || indY > resultMatrix.height()-1)
        {
          break;
        }
    }
  }// end of for loop
}

template<typename T>
void DirectionalLocalExtremSearch( CImg<T>& img, CImg<T>& imgOut, const float threshold, int startX, int startY,
                                  const int gridX, const int gridY, const int dirX, const int dirY)
{
    int Radius = 5;
    int Minimum_Grid = 10;
    int dx = cimg::sign(dirX);
    int dy = cimg::sign(dirY);

    if( (dirX == 0 && dirY == 0) || (gridX < Minimum_Grid) || (gridY < Minimum_Grid) )
    {
        return;
    }

    int stepX = dirX*gridX;
    int stepY = dirY*gridY;

    int centerX = startX;
    int centerY = startY;

    int localMaxX, localMaxY;

    while(centerX>Radius && centerX < img.width()-Radius &&
          centerY>Radius && centerY < img.height()-Radius)
    {
        int count = SearchLocalExtrem(img, centerX, centerY, Radius, localMaxX, localMaxY, threshold);

        if(count > 0)
        {
          int x = localMaxX;
          int y = localMaxY;
          float sum3x3 = img(x-1, y-1) + img(x, y-1) + img(x+1, y-1) +
                         img(x-1, y  ) + img(x, y  ) + img(x+1, y  ) +
                         img(x-1, y+1) + img(x, y+1) + img(x+1, y+1);

          float avgX =  (x*sum3x3  - (img(x-1, y-1) + img(x-1, y)+img(x-1, y+1))
                                   + (img(x+1, y-1) + img(x+1, y)+img(x+1, y+1))
                        )/sum3x3;

           float avgY = (y*sum3x3  - (img(x-1, y-1) + img(x, y-1)+img(x+1, y-1))
                                   + (img(x-1, y+1) + img(x, y+1)+img(x+1, y+1))
                        )/sum3x3;

            imgOut(localMaxX, localMaxY, 0, 0) = 3000;
            imgOut(localMaxX, localMaxY, 0, 1) = avgX;
            imgOut(localMaxX, localMaxY, 0, 2) = avgY;
           
            centerX = localMaxX + stepX;
            centerY = localMaxY + stepY;
        }
        else
        {
            centerX += stepX;
            centerY += stepY;
        }
    }
}


template<typename T>
int
Get_Local_Max(CImg<T>& img, int x0, int y0, int radius,
                            int x1, int y1, int& lx, int& ly, float threshold)
{
  int R = radius;
  CImg<T> output(img);
  output.fill(0);
  int count = 0;
  bool unqualified = false;

  if( R < 0  || x0 > x1 || y0 > y1  ||
      x0 < R || x1 >= img.width() - R ||
      y0 < R || y1 >= img.height() - R )
  {
    return count;
  }

  CImg<T> cropImg = img.get_crop(x0, y0, 0, 0, x1, y1, 0, 0);

  T loc_min, loc_max;
  loc_min = cropImg.min_max(loc_max);

  long maxIndex = 0;
  loc_min = cropImg.min_max_index(loc_max, maxIndex);
 
  cimg_for_insideXY(cropImg, x, y, 0)
  {
    if(cropImg(x,y) == loc_max)
    {
      if( loc_max < threshold )
      {
        unqualified = true;
      }
      else
      {
        lx = x+x0; ly = y+y0;
        count++;
        return count;
      }
    }
  }

  if(count==0)
  {
    if(unqualified)
      count = -1;
  }

  return count;
}

template<typename T>
void
OutputResult(CImg<T>& img)
{
	FILE* pf = fopen("ResultOutput.txt", "w");

	cimg_forY(img, y){
		cimg_forX(img, x)
		{
			if(img(x,y,0,0))
				fprintf(pf,"%7.3f %7.3f\n", img(x, y, 0, 1), img(x, y, 0, 2));
		}
	}
	fprintf(pf, "\n");
	fclose(pf);
}

bool HSCentroid(std::string fileName )
{ 
  CImg<float>  imgIn(fileName.c_str());

  imgIn.blur(4.0);

  int stepX, stepY;
  int cX, cY; // maximum pixel location near center;
  bool bOk = FindSteSize( imgIn, cX, cY, stepX, stepY);

  int startIndexX = cimg::max<float>(cX, imgIn.width()-cX)/stepX + 1;
  int startIndexY = cimg::max<float>(cY, imgIn.height()-cY)/stepY + 1;

  int matrixX = 2*startIndexX + 1;
  int matrixY = 2*startIndexY + 1;

  CImg<float> resMatrix(matrixX, matrixY, 1, 3);
  resMatrix.fill(0);

  if(bOk)
  {
    float global_threshold = 6;
    CImg<float> imgOut(imgIn.width(), imgIn.height(), 1, 3);
    imgOut.fill(0);

    long t0 = cimg::time();

   // DirectionalLocalExtremSearch(imgIn, imgOut, global_threshold, cX, cY, stepX, stepY,  1,  1);

    //DirectionalLocalExtremSearchX( imgIn, resMatrix, global_threshold,
    //                                cX, cY,
    //                                startIndexX, startIndexY,
    //                                stepX, stepY);

	DirectionalLocalExtremSearchY( imgIn, resMatrix, global_threshold,
                                    cX, cY,
                                    startIndexX, startIndexY,
                                    stepX, stepY);

	
    long t1 = cimg::time();

	for(int i = 0; i<resMatrix.height(); i++)
	{
		if(resMatrix(startIndexX, i)>0)
		{
			int cx = (int)resMatrix(startIndexX, i, 0, 1);
			int cy = (int)resMatrix(startIndexX, i, 0, 2);

			DirectionalLocalExtremSearchX( imgIn, resMatrix, global_threshold,
                                    cx, cy,
                                    startIndexX, i,
                                    stepX, stepY);
		}
	}

    long t2 = cimg::time();

    long dt1 = t1 - t0;
    long dt2 = t2 - t1;

	resMatrix.save_bmp("ResultImage.bmp");
	OutputResult(resMatrix);

    resMatrix.display("Result Matrix", true);

    //CImg<float> norm(imgIn.get_structure_tensors());
    //norm.get_normalize(0, 255);
    //norm.display("Gradient", true);

  }

  return bOk;
}

void CImage_Test0(std::string fileName)
{
  CImg<unsigned char>  imgIn(fileName.c_str());
  double mean = imgIn.mean();
  imgIn.blur(2.5);
  mean = imgIn.mean();

  // get buffer data
  const unsigned char* buffer = imgIn.data();

  // create an image instance from given data
  int sizeX = 21, sizeY = 21;
  CImg<float> mask(sizeX, sizeX);

  float dx = 100.0;
  float dy = 100.0;
  float tensor = 4.0;
  float opacity = 0.5;
  const unsigned char white[] = {155, 0, 0};

  int x[2]={0,1200}, y[2] = {1,1000};
  CImg<float> subImage = imgIn.get_crop( x[0], y[0], 0, 0, x[1], y[1], 0, 0);

  float ax=0, ay=0, az=0, ac=0;
  float pixel = imgIn(101,101);

  int x0=200, y0=200, radius=3,
      x1=250, y1=250, lx, ly, count=0;
  float threshold = 30;

 

  count = Get_Local_Max(imgIn, x0, y0, radius, x1, y1, lx, ly, threshold);

  int stepX, stepY;
  int maxX, maxY;

  bool found = FindSteSize ( imgIn, maxX, maxY, stepX, stepY);

  imgIn.display("input image", true);

 // return;

  CImg<float>  subImage2(subImage);
  subImage2.pow(2);

  mask.fill(1);
  subImage2.convolve(mask);

  CImg<float> gaussianMask(mask);
  gaussianMask.fill(0);
  gaussianMask.draw_gaussian(10, 10, tensor, white, opacity);

  CImg<float> gaussianMask2(gaussianMask);
  float Gaussian_Mask_Norm = gaussianMask2.sum();

  
  float flmax = 0;
  float flmin = gaussianMask.min_max(flmax);
  gaussianMask.display("gaussian mask", true);

  CImg<float> convImage(subImage); 
  convImage =  convImage.convolve(gaussianMask);

  flmin = convImage.min_max(flmax);

  subImage2.display("SubImage Pow 2", true);

  convImage.display("2D Gaussian Mask Image before Normalize", true);

  cimg_for_insideXY(convImage,x,y,1) { convImage(x,y) = convImage(x,y)/(sqrt(subImage2(x,y)) + Gaussian_Mask_Norm);}

  flmin = convImage.min_max(flmax);

  convImage.normalize(0, 155);

  found = true;
  threshold = 20;
  for(int j=10; j<1000; j+=36)
  {
    for(int i=10; i<1200; i+=36)
    {
        found = Get_Local_Max(convImage, i, j, radius, i+40, j+40, lx, ly, threshold);
    }
  }

  convImage.display("normalized convolved Image", true);

  subImage.display("Cropped Image", true);
     
  return;
}