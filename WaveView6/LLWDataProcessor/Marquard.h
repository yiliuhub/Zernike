

/*****************************************************************************
 FILE geometry.h
      define the TYPE of features
******************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "DllExport.h"

#define PIE		3.1415926535897932384626433832795
//				3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647

//ANGLE_E and DIST_E where changed from 1.0e-16 to 1.0e-15 because sometimes
//1.000000 - 1.000000 equals 2.22e-16. DON'T CHANGE
/*1.0e-15 difference is relatively small enough for two angles to consider them equal */ 
#define ANGLE_E 1.0e-15
#define DIST_E 1.0e-15

#define DOUBLE_E 1.0e-10 
#define EPS 1.e-100			// a very small number to avoid divide by zero 

/***************************************************************************
 * 2-D features in 3-D space
 ***************************************************************************/

/*-----------------------------*
 | POINT                       |
 *-----------------------------*/

typedef struct {
        double x; /* cartesian coordinate of X-AXIS */
        double y; /* cartesian coordinate of Y-AXIS */
        double z; /* cartesian coordinate of Z-AXIS */
}  Point, *LPPoint;


/*************************************************************************** 
 * coordinate 
 ***************************************************************************/ 
  
typedef struct {
    Point origin;  /* position relative to the universal coordinate */
    Point x_axis;  /* direction of x_axis */
    Point y_axis;  /* direction of y_axis */
    Point z_axis;  /* direction of z_axis */
                       /* all directions relative to the universal coordinate */ 
} Coordinate;

/*-----------------------------*
 | LINE                        |
 *-----------------------------*/

typedef struct {
    Point  middle_point; /* middle point of the line */
    double length;       /* length of the line */
    Point  direction;    /* unit direction vector which lies along the line 
                            and points from the first point to the second */
} Line; 

/*-----------------------------*
 | CIRCLE                      |
 *-----------------------------*/

typedef struct {
    Point  center;    /* center of the circle */
    Point  direction; /* unit direction vector of the plane where circle lies in*/
    double radius;	  /* radius of the circle */
} Circle;




/*-----------------------------------------------* 
 | HYPERBOLA --- Intersection of CONE and PLANE  | 
 *-----------------------------------------------*/ 
 
typedef struct {
    Point  first_focus;  /* one focus of the hyperbola */ 
    Point  second_focus; /* one focus of the hyperbola */ 
    double radius;       /* major radius in which two focuses lie */ 
    Point  direction;    /* unit direction vector of the plane where Hyperbola lies in*/
} Hyperbola; 


/*-----------------------------------------------*  
 | PARABOLA --- Intersection of CONE and PLANE   |  
 *-----------------------------------------------*/  
  
typedef struct { 
    Point  focus;     /* one focus of the Parabola */  
    Line   directrix; /* directrix of the Parabola */  
    Point  direction; /* unit direction vector of the plane where the Parabola lies in*/
} Parabola;  

/*-----------------------------*
 | PLANE                       |
 *-----------------------------*/
 
typedef struct {
        Point  center;    /* a point on the plane */
        Point  direction; /* unit direction vector of the plane */ 
} Plane; 
 

/*-----------------------------*
 | CONE                        |
 *-----------------------------*/

typedef struct {
    Point  vertex;       /* vertex of the cone */
    double angle;        /* included angle of the cone */
    double first_height; /* distances from the vertex to two parallel planes which */
    double last_height;  /* perpendiculer to corn's direction and cut the cone*/
    Point  direction;    /* unit direction vector associated with cone, which points 
                            along the cone's axis from the vertex to the open end */
} Cone;

/*-----------------------------*
 | CYLINDER                    |
 *-----------------------------*/
 
typedef struct {
    Point  center;       /* center of the cylinder */
    double radius;       /* radius of the cylinder */ 
    Point  direction;    /* unit direction vector associated with the cylinder, and
                            points along the cylind's axis from the first end to the 
                            last end */
} Cylinder;

/*-----------------------------*
 | SPHERE                      |
 *-----------------------------*/
  
typedef struct {
    Point  center;   /* center of the Sphere */
    double radius;   /* radius of the Sphere */ 
} Sphere;
 
/*----------------------------------------------------------------*
 | GAUSSIAN FUNCTION :											  |
 |				amp * exp( -(x-ctr)*(x-ctr)/(dev*dev)			  |
 |																  |
 *----------------------------------------------------------------*/
  
typedef struct {
    double amp;		//amplifier
    double ctr;		//center
	double dev;		//deviation
} Gaussian;

/*----------------------------------------------------------------*
 | 2D GAUSSIAN FUNCTION 
 | Type 1:														|
 |			amp * exp(	-(x-ctrX)*(x-ctrX)/(devX*devX) 
 |						-(y-ctrY)*(y-ctrY)/(devY*devY)
 |					 )											  |
 |																  |
 *----------------------------------------------------------------*/
typedef struct {
    double amp;			//amplifier
    double ctrX;		//center
	double ctrY;		//center
	double devX;		//deviation
	double devY;		//deviation
} Gaussian2D;


/*----------------------------------------------------------------------------------*
 | 2D GAUSSIAN FUNCTION																|
 | Type 2:
 | G(x, y) = amp * exp( a*(x-ctrX)^2 + b*(x-ctrX)*(y-ctrY) + c*(y-ctrY)^2 )			|
 |																					|
 *----------------------------------------------------------------------------------*/

typedef struct {
    double amp;		//amplifier
    double ctrX;	//center
	double ctrY;	//center
	double a;		//coefficient of x^2
	double b;		//coefficient of x*y
	double c;		//coefficient of y^2
} GaussianAffine2D;

/*-------------------------------------------------------------------------*
 |  Marquard Numerical Optimization
 *-------------------------------------------------------------------------*/
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
   double   *sum_sq);

/*-------------------------------------------------------------------------*
 |  Data performance
 *-------------------------------------------------------------------------*/
Point v_create(double x, double y, double z);
double distance_point_point(Point p1, Point p2);

/*-------------------------------------------------------------------------*
 |  Newton Minimization 
 *-------------------------------------------------------------------------*/
int newton_minimize(double (*fun_der1)(double*, double), double (*fun_der2)(double*, double),
							  double *parameters, double *minimizer, double init_value, double tolerence);

/*---------------------------------------------------------------*
 | GAUSSIAN --- implementation is in gaussian.c                  |
 *---------------------------------------------------------------*/
int gaussian_minimizer_point(double* point, double* b, Point* pt, double* x);
double gaussian_sum_fun(double* point, double* b);
double point_to_gaussian_sq_D1(double *param, double x);
double point_to_gaussian_sq_D2(double *param, double x);
int initial_gaussian(Point* point, int n, Gaussian* init_gaussian);
int gaussian_bestfit(Point *points, int n, double tolerance, int iterateLimit, 
					Gaussian *init_gaussian, int useInit, double *extraInfo);
double gaussian(double amf, double ctr, double dev, double x);


int gaussian2D_minimizer_point(double* point, double* b, Point* pt, double* x);
double gaussian2D_sum_fun(double* point, double* b);
double point_to_gaussian2D_sq_D1(double *param, double x);
double point_to_gaussian2D_sq_D2(double *param, double x);
int initial_gaussian2D(Point* point, int n, Gaussian* init_gaussian);
int gaussian2D_bestfit(Point *points, int n, double tolerance, int iterateLimit, 
						Gaussian2D *init_gaussian, int useInit, double *extraInfo);
double gaussian2D(double amf, double ctrX, double ctrY, double devX, double devY, double x, double y);


int gaussianAffine2D_minimizer_point(double* point, double* b, Point* pt, double* x);
double gaussianAffine2D_sum_fun(double* point, double* b);
double point_to_gaussianAffine2D_sq_D1(double *param, double x);
double point_to_gaussianAffine2D_sq_D2(double *param, double x);
int initial_gaussianAffine2D(Point* point, int n, GaussianAffine2D* init_gaussian);
int gaussianAffine2D_bestfit(Point *points, int n, double tolerance, int iterateLimit, 
						GaussianAffine2D *init_gaussian, int useInit, double *extraInfo);
double gaussianAffine2D(double amf, double ctrX, double ctrY, double a, double b, double c, double x, double y);

double gaussianAffine2D_deriv_fun(double* point, double* b, double* d);


#ifdef __cplusplus
}
#endif

#endif

