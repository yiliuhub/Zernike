
/**********************************************************************
 * File: Gaussian.c
 *
 * Yi Liu 01/10/02
 **********************************************************************/

#include <math.h>

#include "Geometry.h"
#include "geoerror.h"

double gaussian(double amp, double ctr, double dev, double x)
{
	double value = 0;

	if(dev*dev <= DOUBLE_E)	// the dev is 0
		value = 0;			// return 0 to punish the invalid deviation
	else
	{
        value = amp*exp(-(x - ctr)*(x - ctr)/(dev*dev));
	}

	return value;
}


/**********************************************************************
*	The parameterized equation of a gaussian is:
*		P(x) = A * exp( -(x-C)^2/ (B*B))
*
*	here 
*		A --- amplifier
*		C --- center
*		B --- deviation
*
***********************************************************************/

int gaussian_minimizer_point(double* point, double* b, Point* pt, double* x)
{
	int result;
	double A = 1, B = 1, C = 0;
	double param[5];
	double start_value;
	double num_tolerence = 1.e-20;
	
	start_value = point[0];
	param[0] = point[0];
	param[1] = point[1];
	param[2+0] = b[0];
	param[2+1] = b[1];
	param[2+2] = b[2];
	
	result = newton_minimize(point_to_gaussian_sq_D1, point_to_gaussian_sq_D2,
								param, x, start_value, num_tolerence);

	A = param[2];
	B = param[3];
	C = param[4];

	// the minimizing point on ellipse
	pt->x = *x;
	pt->y = A*exp(-(*x - C)*(*x - C)/(B*B));
	pt->z = 0;

	return result;
}

double gaussian_sum_fun(double* point, double* b)
{
	int result;
	double s;
	Point pt;
	Point p;

	p = v_create(point[0], point[1], 0);
	result = gaussian_minimizer_point(point, b, &pt, &s);

	//gaussian is a one variable function
	pt.z = 0;

	return distance_point_point(p, pt);
}


/******************************************************************************************
*
Gaussian like function of variable x is expressed as

gauss(x) = A*exp(-(x-C)*(x-C)/(B*B))

here A, B, and C are parameters. The fitting of a set of points with Gaussian like function \
needs to calculate the partial derivtives with respect to A, B, and C.

*  Yi Liu 01/16/01
*
******************************************************************************************/
static double gaussian_deriv_fun(double* point, double* b, double* d)
{
	double x, x0, y0;
	double A, B, C;
	double D;
	double gauss;
	Point Pt;
	

	gaussian_minimizer_point(point, b, &Pt, &x);
	A = b[0];
	B = b[1];
	C = b[2];
	x0 = point[0];
	y0 = point[1];

	gauss = A*exp(-(x-C)*(x-C)/(B*B));
	D = (x-x0)*(x-x0) + (gauss-y0)*(gauss-y0);
	D = sqrt(D);

	if(D > 0)
	{
		//partial derivative with respect to A (or amplifier)
		d[0] = gauss*(gauss - y0)/(A*D);

		//partial derivative with respect to B (or deviation)
		d[1] = -2*gauss*(gauss-y0)*(x-C)*(x-C)/(B*B*B*D);

		//partial derivative with respect to C (or center)
		d[2] = 2*(x-C)*gauss*(gauss-y0)/(B*B*D);
	}
	else
	{
		d[0] = 0;
		d[1] = 0;
		d[2] = 0;
	}

	return 1.0;
}


int initial_gaussian(Point* points, int n, Gaussian* init_gaussian)
{
	int i;
	int count;
	double DevSQ = 0;
	double C;
	double A;
	
	A = points[0].y;
	C = points[0].x;

	for(i = 1; i < n; i++)
	{
		if(fabs(points[i].y) > A)
		{	
			A = points[i].y;
			C = points[i].x;
		}
	}

	DevSQ = 0;
	count = n;
	for(i = 0; i < n; i++)
	{
		if(points[i].y/A < 1)
			DevSQ += - (points[i].x - C)*(points[i].x - C)/(log(points[i].y/A));
		else
			count--;
	}

	DevSQ /= count;

	init_gaussian->amp = A;
	init_gaussian->ctr = C;
	init_gaussian->dev = sqrt(DevSQ);

	return 1;
}


int gaussian_bestfit(Point *points, int n, double tolerance,
					int iterateLimit, Gaussian *init_gaussian,
					int useInit, double *extraInfo)
{
   double solVec[5];
   int success;
 
   if (n < 4)
   {
      error_coding(TOO_FEW_POINTS);
      return FAILURE;
   }


   if(!useInit)
   {
		success = initial_gaussian(points, n, init_gaussian);

		if (!success)
			return FAILURE;
   }

   
   solVec[0] = init_gaussian->amp;
   solVec[1] = init_gaussian->dev;
   solVec[2] = init_gaussian->ctr;

   success = marquardt((double*)points, n, gaussian_sum_fun, gaussian_deriv_fun, 
                       solVec, 3, tolerance, iterateLimit, extraInfo);

   if (!success)
      return FAILURE;

   init_gaussian->amp = solVec[0];
   init_gaussian->dev = solVec[1];
   init_gaussian->ctr = solVec[2];
   
   return success;
}

/****************************************************************
* point_to_gaussian_sq_D1() --- the first derivative of the distance 
*							   function from point to an point on 
*							   gaussian curv
*****************************************************************/
double point_to_gaussian_sq_D1(double *param, double x)
{
	double A, B, C;
	double B2;
	double x0, y0;
	double gauss;
	double derivative;

	x0 = param[0];
	y0 = param[1];

	A = param[2];
	B = param[3];
	C = param[4];

	B2 = B*B;
	gauss = A*exp(-(x-C)*(x-C)/B2);

	derivative = 2*(x - x0) - 4*(x-C)/B2 * gauss*(gauss - y0);

	return derivative;
}



/******************************************************************
* point_to_gaussian_sq_D2() --- the second derivative of the distance 
*							   function from point to a point on
*								gaussian curv
*
*******************************************************************/
double point_to_gaussian_sq_D2(double *param, double x)
{
	double A, B, C;
	double B2;
	double x0, y0;
	double gauss;
	double derivative;

	x0 = param[0];
	y0 = param[1];

	A = param[2];
	B = param[3];
	C = param[4];

	B2 = B*B;
	gauss = A*exp(-(x-C)*(x-C)/B2);

	derivative = 2 + y0*(2/B2 - 4*(x-C)*(x-C)/(B2*B2))*gauss
				   + (-4/B2 + 16*(x-C)*(x-C)/(B2*B2))*gauss*gauss;
					

	return derivative;
}


/*****************************************************************************************
*
*	2D Gaussian Bestfit functions
*
*	A 2D normalized Gaussian function looks like
*						1
*	G(x,y;a,b) =  ----------------- exp(-(x-x0)^2/(2*a) - (y-y0)^2/(2*b))
*					2*PI* sqrt(a*b)
*
*	If we want a non-nomralized 2D Gaussian function, it will look like
*
*
*	(*)		G(x,y;a,b) =  C * exp(-(x-x0)^2/(2*A) - (y-y0)^2/(2*B))
*
*	In real cases, data to be fitted are not from normalized (Gaussian) functions, so 
*	we'd adapt to (*)
*
****************************************************************************************/






/*****************************************************************************
*	Name: gaussian2D_minimizer_point
*	Function: Calculate the distance from *point to a 2D Gaussian function
*
*	The parameterized 2D gaussian function is:
*
*		G(x,y;A,B,C) =  C * exp(-(x-x0)^2/(2*A) - (y-y0)^2/(2*B))
*
*	here 
*		(x0, y0) --- center
*		C --- amplifier
*		A --- deviation in x direction
*		B --- deviation in y direction
*
*	The point on the surface of 2D Gaussian function with x, y is:
*		(x, y, G(x, y; A, B, C)
*
*	So what this function does is to find the minimum of the distance 
*	from *point to (x, y, G(x, y; A, B, C)) for all (x, y)
*
*	Yi Liu 01/24/2007
******************************************************************************/

int gaussian2D_minimizer_point(double* point, double* b, Point* pt, double* x)
{
	int result;
	double A = 0, B = 0, C = 1;
	double param[5];
	double start_value;
	double num_tolerence = 1.e-20;
	
	start_value = point[0];
	param[0] = point[0];
	param[1] = point[1];
	param[2+0] = b[0];
	param[2+1] = b[1];
	param[2+2] = b[2];
	
	result = newton_minimize(point_to_gaussian_sq_D1, point_to_gaussian_sq_D2,
								param, x, start_value, num_tolerence);

//	To Be comleted
//	result = steepest_minimize();

	A = param[2];
	B = param[3];
	C = param[4];

	pt->x = *x;
	pt->y = A*exp(-(*x - C)*(*x - C)/(B*B));
	pt->z = 0;

	return result;
}


int steepest_minimize()
{
	return 0;
}

/*
double SteepestOptimization(double* pDist, double* pGradX, double* pGradY, 
				int count, double* pX, double* pY, double& lambda, double& theta, double& Tx, double& Ty)
{
	int i = 0;
	int j = 0;
	int k = 0;
	int iter_num = 1000;
	double min_sum = 1.e60;
	double min_theta = theta, min_Tx = Tx, min_Ty = Ty;
	double sum = 0, new_sum = 0, sumX = 0, sumY =0;
	double gradTheta = 0, gradTx = 0, gradTy = 0;
	double new_theta, new_Tx, new_Ty;
	double t = 10, t0 = 10;

	for(i = 0; i < count; i++)
	{
		Cx += pX[i];
		Cy += pY[i];
	}

	Cx /= count;
	Cy /= count;

	for(i = 0; i < count; i++)
	{
		sum += sqrt((pX[i]-Cx)*(pX[i]-Cx) + (pY[i]-Cy)*(pY[i]-Cy));
	}

	sum /= count;

	if(sum > 1)
		lambda = 2/sum;
	else
		lambda = 1;

	//initialization
	

	for(k = 0; k <= 4; k++)
	{
		if(k < 4)
		{
			theta = min_theta + k*3.14/(2*lambda);
			Tx = min_Tx;
			Ty = min_Ty;
			sum = 1.e60;
		}
		else
		{
			theta = min_theta;
			Tx = min_Tx;
			Ty = min_Ty;
			sum = min_sum;
			t0 = 1;
		}


		for(i = 0; i < iter_num; i++)
		{
			sum = Sum(width, height, pDist, count, pX, pY, lambda, theta, Tx, Ty);

			CalcGradientVector( width, height, pDist, pGradX, pGradY, count,
								pX, pY, lambda, theta, Tx, Ty, gradTheta, gradTx, gradTy);

			new_theta = theta - t*gradTheta;
			new_Tx = Tx - t*gradTx;
			new_Ty = Ty - t*gradTy;

			t = t0;
			for(j = 0; j < 10; j++)
			{
				new_sum = Sum(width, height, pDist, count, pX, pY, lambda, new_theta, new_Tx, new_Ty);

				if(new_sum < sum)
				{
					theta = new_theta;
					Tx = new_Tx;
					Ty = new_Ty;
					break;
				}
				else
				{
					t = 0.3*t;
					new_theta = theta - t*gradTheta;
					new_Tx = Tx - t*gradTx;
					new_Ty = Ty - t*gradTy;
				}
			}

	 		if(j == 10)
				break;
		}

		if(min_sum > sum)
		{
			min_sum = sum;
			min_theta = theta;
			min_Tx = Tx;
			min_Ty = Ty;
		}
	}
	
	Tx = min_Tx;
	Ty = min_Ty;
	theta = min_theta;

	return min_sum;
}
*/
