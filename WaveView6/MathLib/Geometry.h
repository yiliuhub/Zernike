

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

/*-----------------------------*
 | ARC                         |
 *-----------------------------*/

typedef struct {
    Point  center;      /* center of the circle where the arc lies along */
    double radius;      /* radius of the circle where the arc lies along */
    Point  direction;   /* unit direction vector of the plane where arc lies in*/ 
    Point  start_point; /* start point of the arc */
    Point  end_point;   /* end point of the arc; by right-hand rule */
} Arc;

/*---------------------------------------------*
 | ELLIPSE --- Intersection of CONE and PLANE  |
 *---------------------------------------------*/

typedef struct {
    Point  first_focus;  /* one focus of the ellipse */
    Point  second_focus; /* one focus of the ellipse */
    double radius;       /* major radius in which two focuses lie */
    Point  direction;    /* unit direction vector of the plane where ellipse lies in*/
} Ellipse;

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
 | WIDTH                       |
 *-----------------------------*/
  
typedef struct {
    Point  location;	/* midpoint of the width	*/
    double width;		/* the actual width			*/ 
	double angle;		/* center line angle		*/
} WidthStruc;
 
/*-----------------------------*
 | CONTOUR                      |
 *-----------------------------*/

typedef struct {
    Point  center;    /* center of the contour */
    Point  direction; /* unit direction vector of the plane that contour lies in*/
    double length;	  /* length of the contour */
	double area;	  /* area of the contour   */
} Contour;

/***************************************************************************
 * 3-D features in 3-D space
 ***************************************************************************/

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


/********************************************************************************
 * Global Variables 
********************************************************************************/

extern DllExport Point V_I;
extern DllExport Point V_J;
extern DllExport Point V_K;
extern DllExport Point V_0;

/********************************************************************************
 * Prototype of functions
********************************************************************************/

/*---------------------------------------------------------------*
 | MEMORY --- implemention is in memory.c                        |
 *---------------------------------------------------------------*/
DllExport short* new_space_short( int n );
DllExport int* new_space_int( int n );
DllExport float* new_space_float( int n );
DllExport double* new_space_double( int n );
DllExport Point* new_space_point( int n );
DllExport Line* new_space_line( int n );
DllExport Circle* new_space_circle( int n );
DllExport Arc* new_space_arc( int n );
DllExport Ellipse* new_space_ellipse( int n );
DllExport Parabola* new_space_parabola( int n );
DllExport Hyperbola* new_space_hyperbola( int n );
DllExport Cone* new_space_cone( int n );
DllExport Sphere* new_sphere( int n );
DllExport Cylinder* new_space_cylinder( int n );

/*---------------------------------------------------------------*
 | DISTANCE --- implemention is in distance.c                    |
 *---------------------------------------------------------------*/
DllExport double distsq_circle_point(Circle circle, Point point) ;
DllExport double distsq_line_line( Line line1, Line line2  );
DllExport double distsq_line_plane( Line line, Plane plane );
DllExport double distsq_line_point( Line line, Point point );
DllExport double distsq_plane_plane( Plane plane1, Plane plane2 );
DllExport double distsq_plane_point( Plane plane, Point point );
DllExport double distsq_point_point( Point point1, Point point2);
DllExport double distsq_point_sphere( Point point, Sphere sphere );
DllExport double distsq_cone_point( Cone cone, Point point );
DllExport double distsq_cylinder_point( Cylinder cylinder, Point point );

DllExport double distance_circle_point(Circle circle, Point point );
DllExport double distance_line_line( Line line1, Line line2 );
DllExport double distance_line_plane( Line line, Plane plane );
DllExport double distance_line_point( Line line, Point point );
DllExport double distance_plane_plane( Plane plane1, Plane plane2 );
DllExport double distance_plane_point( Plane plane, Point point );
DllExport double distance_point_point( Point point1, Point point2);
DllExport double distance_point_sphere( Point point, Sphere sphere );
DllExport double distance_cone_point( Cone cone, Point point );
DllExport double distance_cylinder_point( Cylinder cylinder, Point point );

/*---------------------------------------------------------------*
 | RELATION --- implementation is in relation.c                  |
 *---------------------------------------------------------------*/
DllExport int parallel_directions( Point direction1, Point direction2 );
DllExport int parallel_line_line( Line line1, Line line2 );
DllExport int parallel_line_plane( Line line, Plane plane );
DllExport int parallel_plane_plane( Plane plane1, Plane plane2 );

DllExport int perpendicular_directions( Point direction1, Point direction2 );
DllExport int perpendicular_line_line( Line line1, Line line2 );
DllExport int perpendicular_line_plane( Line line, Plane plane );
DllExport int perpendicular_plane_plane( Plane plane1, Plane plane2 );

DllExport int overlap_line_line( Line line1, Line line2 );
DllExport int overlap_line_plane( Line line, Plane plane );
DllExport int overlap_plane_plane( Plane plane1, Plane plane2 );
DllExport int overlap_point_line( Point point, Line line );
DllExport int overlap_point_plane( Point point, Plane plane );

DllExport int convex_points(Point* points, int n, Point* workspace);

/*---------------------------------------------------------------*
 | ANGLE_BETWEEN --- implementation is in angle_between.c        |
 *---------------------------------------------------------------*/
DllExport double cos_angle_directions( Point direction1, Point direction2 ) ;
DllExport double angle_directions( Point direction1, Point direction2, 
                         Point* cross_direction );

/*---------------------------------------------------------------*
 | MATRIX --- implementation is in matrix.c                      |
 *---------------------------------------------------------------*/
DllExport double* get_cofactor_matrix( double* matrix, int n, int i, int j, 
                             double* cofactor_matrix);
DllExport double determinant( double* matrix, int n );
DllExport double* cramer_rule_matrix( double* matrix, double* vector, int n, int j, 
                            double* cramer_matrix );
DllExport double* liner_equation( double* matrix, double* vector, int n, 
                        double* solution );

 /*---------------------------------------------------------------*
 | MATH_OPERATION --- implementation is in math_operation.c      |
 *---------------------------------------------------------------*/
DllExport Point v_add( Point vector1, Point vector2 );
DllExport Point v_sub( Point vector1, Point vector2 );
DllExport double v_dot( Point vector1, Point vector2 );
DllExport Point v_mult_c( Point vector, double coefficient );
DllExport Point v_cross( Point vector1, Point vector2 );
DllExport Point* v_copy( Point* v_addr1, Point* v_addr2 );
DllExport int v_equal( Point vector1,  Point vector2 );
DllExport int v_equal_e( Point vector1,  Point vector2, double e );
DllExport Point v_create( double x, double y,  double z );
DllExport double v_distsq( Point point );
DllExport double v_dist( Point point );
DllExport Point v_center( Point point1, Point point2  );
DllExport int d_equal( double x, double y, double e );
DllExport int quadratic_equation( double a, double b, double c, 
                        double* root1, double* root2 );

DllExport int merge_sort(double* values, int n, int* indexes);

/*---------------------------------------------------------------*
 | COORDINATE --- implementation is in coordinate.c              |
 *---------------------------------------------------------------*/
DllExport Point point_coordinate_move(Coordinate old_coord, Point point, 
								Coordinate new_coord);
DllExport Point point_coordinate_transform(Coordinate old_coord, Point point, 
                                Coordinate new_coord);
DllExport Point unit_direction(Point direction);
DllExport void normal_coordinate(Coordinate* coordinate);
DllExport Coordinate coordinate_direction (Coordinate coordinate, Point direction, 
                                                        int axis);
DllExport void coordinate_transform(Coordinate old_coord, double* points, int n, 
								Coordinate new_coord);
DllExport Line line_coordinate_transform(Coordinate new_coord, Line line, 
								Coordinate old_coord);
DllExport Circle circle_coordinate_transform(Coordinate old_coord, Circle circle, 
                                Coordinate new_coord) ;

/*---------------------------------------------------------------*
 | PROJECTION --- implementation is in projection.c              |
 *---------------------------------------------------------------*/
DllExport Point point_direction_move( Point point, Point direction, double distance );
DllExport Point projection_point_line( Point point, Line line );
DllExport Point projection_point_circle( Point point, Circle circle );
DllExport Point projection_point_plane( Point point, Plane plane );
DllExport Line  projection_line_plane( Line line, Plane plane );
DllExport int projection_points_plane(Point* points, int n, Plane plane ) ;

/*---------------------------------------------------------------*
 | INTERSECTION --- implementation is in intersection.c          |
 *---------------------------------------------------------------*/
DllExport int intersection_line_line( Line line1, Line line2, Point *d_point1, 
                      Point *d_point2, double *angle, Point* cross_direction );

DllExport int intersection_circle_line(Circle circle, Line line, 
                      Point* p_intersect1, Point* p_intersect2 );

DllExport int new_intersection_circle_line(Circle circle, Line line, 
                      Point* p_intersect1, Point* p_intersect2 );

DllExport int intersection_arc_line(Arc arc, Line line, 
                      Point* p_intersect1, Point* p_intersect2 );

DllExport int intersection_arc_line_2D(Arc arc, Line line, 
                      Point* p_intersect1, Point* p_intersect2 );

DllExport int intersection_line_plane( Line line, Plane plane, Point* intersect_addr, 
                double* angle, double* normal_angle, Point* cross_direction );

DllExport int intersection_cone_line( Cone cone, Line line, 
                    Point* p_intersect1, Point* p_intersect2, int* tangent);

DllExport int interLines_cone_plane( Cone cone, Plane plane, Line* intersectLine1, 
                                Line*  intersectLine2, double* include_angle );

DllExport int interEllipse_cone_plane( Cone cone, Plane plane, Ellipse *intersectEllipse, 
                      Circle *intersectcircle, Point *center, double *b_axis );

DllExport int interParabola_cone_plane( Cone cone, Plane plane, 
                               Parabola *intersectParabola, Point *vertex );

DllExport int interHyperbola_cone_plane( Cone cone, Plane plane, 
              Hyperbola* intersectHyperbola, Point* vertex1, Point* vertex2 );

DllExport int intersection_cone_plane(Cone cone, Plane plane, void **intersect_addr );
 
DllExport int intersection_line_sphere( Line line, Sphere sphere, 
                            Point* p_intersect1, Point* p_intersect2 );

DllExport int intersection_plane_sphere( Plane plane, Sphere sphere, 
                             Circle* interCircle, Point* interPoint );

DllExport int intersection_plane_plane( Plane plane1, Plane plane2, Line* interLine, 
                                double* angle, double* normal_angle );

DllExport int intersection_cylinder_line( Cylinder cylinder, Line line,
                                Point *p_intersect1, Point *p_intersect2 );

DllExport int intersection_cylinder_plane( Cylinder cylinder, Plane plane, 
                       Ellipse *interEllipse, Circle* interCircle, 
                       Line *interLine1, Line *interLine2 );

DllExport int intersection_circle_plane(Circle circle, Plane  plane,
							Point *p_intersect1, Point  *p_intersect2);

DllExport int intersection_arc_plane(Arc arc, Plane  plane,
							Point *p_intersect1, Point  *p_intersect2);

/*---------------------------------------------------------------*
 | GAUSSIAN --- implementation is in gaussian.c                  |
 *---------------------------------------------------------------*/
DllExport int gaussian_minimizer_point(double* point, double* b, Point* pt, double* x);
DllExport double gaussian_sum_fun(double* point, double* b);
DllExport double point_to_gaussian_sq_D1(double *param, double x);
DllExport double point_to_gaussian_sq_D2(double *param, double x);
DllExport int initial_gaussian(Point* point, int n, Gaussian* init_gaussian);
DllExport int gaussian_bestfit(Point *points, int n, double tolerance, int iterateLimit, 
							   Gaussian *init_gaussian, int useInit, double *extraInfo);
double gaussian(double amf, double ctr, double dev, double x);


/*---------------------------------------------------------------*
 | POINT --- implementation is in point.c                        |
 *---------------------------------------------------------------*/
DllExport int  check_zero_direction(Point direction);
DllExport int average_dir(Point* points, int n, double* dir);
DllExport double* average_points_2D(Point* points, int n, double* average_p );
DllExport double* average_points(Point* points, int n, double* average_p );
DllExport int average_dir_2D(Point* points, int n, double* dir );
DllExport Point* nearest_point_2_line(Point* points, int n, Line line, Point* near_points );
DllExport Point* nearest_point(Point* points, int n, Point point, Point* near_points );
DllExport double center_moment(Point *points, int n);
 
/*---------------------------------------------------------------*
 | LINE --- implementation is in line.c                          |
 *---------------------------------------------------------------*/
DllExport void end_points_line(Line line, Point* p_point1, Point* p_point2 );
DllExport int line_2_points(Point point1, Point point2, Line* line );

DllExport double farpoint_line_2D(double* points, int n, double* solVec, 
                         double* farpoints );
DllExport double farpoint_line(double* points, int n, double* solVec, 
                         double* farpoints );

DllExport int line_bestfit_2D_calculated(Point* points, int n, Line* init_line, 
                           double* extraInfo);
DllExport int line_bestfit_2D(Point* points, int n, double tolerance, int iterateLimit, 
                          Line* init_line, int useInit, double* extraInfo);
DllExport int line_max_min(Line line, Point* points, int n,
             Line* max_min_line, int max, int*index);

DllExport void line_maximum(Line line, Point* points, int n,
               Line* maxLine, int* index);
DllExport void line_minimum(Line line, Point* points, int n,
               Line* minLine, int* index);

DllExport int line_geom_max_min(Line line, Point* points, int n,
             Line* max_min_line, int max, int*index);

DllExport void line_geometry_maximum(Line line, Point* points, int n,
               Line* maxLine, int* index);

DllExport void line_geometry_minimum(Line line, Point* points, int n,
               Line* minLine, int* index);

DllExport int line_bestfit(Point* points, int n, double tolerance, int iterateLimit, 
                          Line* init_line, int useInit, double* extraInfo);
DllExport int line_fit_within_plane(Point* points, int n, Point direction, 
                                   Line* init_line, double* extraInfo);
DllExport int collinear_points( Point* points, int n);
DllExport int collinear_points_2D( Point* points, int n);

/*---------------------------------------------------------------*
 | ARC --- implementation is in arc.c                            |
 *---------------------------------------------------------------*/
DllExport int arc_3_points(Point point1, Point point2, Point point3, Arc* p_arc);
DllExport double  included_angle_arc(Arc arc );
DllExport int arc_bestfit( Point* points, int n, double tolerance, int iterateLimit, 
                         Arc* init_arc, int useInit, double* extraInfo);
DllExport int arc_bestfit_2D(Point* points, int n, double tolerance, int iterateLimit, 
                        Arc* init_arc, int useInit , double* extraInfo);
DllExport int arc_fit_within_plane( Point* points, int n, Point direction, 
                   double tolerance, int iterateLimit, 
				   Arc* init_arc, int useInit, double* extraInfo);
DllExport int arc_project_bestfit(Point* points, int n, double tolerance, 
              int iterateLimit, Arc* init_arc, int useInit, double* extraInfo);

/*---------------------------------------------------------------*
 | CIRCLE --- implementation is in circle.c                      |
 *---------------------------------------------------------------*/
DllExport double circle_sum_fun(double* point, double* b );
 
DllExport int  average_circle_2D(Point* points, int n, double* average_circle);
DllExport int  average_circle(Point* points, int n, double* average_circle);

DllExport int  initial_circle(Point* point, int n, Circle* raw_circle);

DllExport int circle_3_points( Point point1, Point point2, Point point3, Circle* circle);
DllExport int circle_bestfit( Point* points, int n, double tolerance, int iterateLimit, 
                       Circle* init_circle, int useInit, double* extraInfo);
DllExport int circle_bestfit_2D(Point* points, int n, double tolerance, int iterateLimit, 
                       Circle* init_circle, int useInit, double* extraInfo );
DllExport int circle_fit_within_plane( Point* points, int n, Point direction, 
                 double tolerance, int iterateLimit, 
				 Circle* init_circle, int useInit, double* extraInfo);
DllExport int circle_project_bestfit(Point* points, int n, double tolerance, 
                 int iterateLimit, Circle* init_circle, 
				 int useInit, double* extraInfo);

DllExport void circle_max_min(Circle circle, Point* points, int n, 
               Circle* maxCircle, int maximum, int*index);
DllExport void circle_maximum(Circle circle, Point* points, int n, 
               Circle* maxCircle, int* index);
DllExport void circle_minimum(Circle circle, Point* points, int n, 
               Circle* minCircle, int* index);
DllExport void new_circle_max_min(Circle circle, Point* points, int n, 
						          Circle* maxMinCircle, int maxMinFlag, int* index);

DllExport int searchArcPoint(Point* points, Point* p, int n, Point direction, 
							 Circle* tempCircle, int maxMinFlag, int *pIndex);
DllExport void expendShrinkCircle(Point *points, Point *p, int n, 
					 Point direction, Circle *circle, int maxMinFlag);

DllExport void maxCircle3Points(Point* p, Circle* circle);

DllExport double nonZero(double numerator);
DllExport int acuteTriangle(Point *p);


/*---------------------------------------------------------------*
 | ELLIPSE --- implementation is in ellipse.c                    |
 *---------------------------------------------------------------*/
DllExport int ellipse_bestfit(Point *points, int n, double tolerance, int iterateLimit, 
							  Ellipse *init_ellipse, int useInit, double *extraInfo);
DllExport int ellipse_minimizer_point(double* point, double* b, Point* pt, double* s);
DllExport double ellipse_sum_fun(double* point, double* b);
DllExport double point_to_ellipse_sq_D1(double *param, double t);
DllExport double point_to_ellipse_sq_D2(double *param, double t);
DllExport int initial_ellipse(Point* point, int n, Ellipse* init_ellipse);



/*---------------------------------------------------------------*
 | PLANE --- implementation is in plane.c                        |
 *---------------------------------------------------------------*/
 DllExport int average_plane(Point* points, int n, double* average_plane);
DllExport int plane_3_points(Point point1, Point point2, Point point3, Plane* plane);
DllExport int plane_line_point(Line line, Point point, Plane* plane);
DllExport int plane_bestfit(Point* points, int n, double tolerance, int iterateLimit, 
                Plane* init_plane, int useInit, double* extraInfo);
DllExport int colplaner_points(Point* points, int n);
DllExport int plane_maximum(Plane plane, Point* points, int n, Plane* maxPlane,
                 double* profile_plus, double* profile_minus, int* indexes);
DllExport int plane_minimum(Plane plane, Point* points, int n, Plane* minPlane,
                 double* profile_plus, double* profile_minus, int* indexes);

DllExport int max_min_plane(Plane plane, Point* points, int n, Plane* maxMinPlane, int maxMinFlag,
							double* profile_plus, double* profile_minus, int* p_indexes);

DllExport int pointInsideTriangle( Point point, Point p1, Point p2, Point p3);
DllExport int searchPoint( Point p1, Point spinAxis, Point normal,
				 Point posiDirect, Point* pts, int n);

DllExport Point* plane_touchProbe_correction(Plane plane, Point* points, int n, 
											 double probeRadius, int getIn);

DllExport int plane_bestfit_touchProbe(Point *points, int n, double tolerance,
         int iterateLimit, Plane  *init_plane, int useInit, double *extraInfo, 
         double probeRadius, int getIn, Point** correctedPoints);

/*---------------------------------------------------------------*
 | SPHERE --- implementation is in sphere.c                      |
 *---------------------------------------------------------------*/
DllExport double sphere_sum_fun(double* point, double* b );

DllExport int average_sphere(Point* points, int n, double* average_sphere );
 
DllExport int sphere_4_points(Point* four_points, Sphere* sphere ); 
DllExport int sphere_bestfit( Point* points, int n, double tolerance, int iterateLimit, 
                    Sphere* init_sphere, int useInit, double* extraInfo);
DllExport int searchSurfacePoint(Point* points, int n, Point sidePoint, Point midPoint,
					   Point direction, Point* objectCenter, int maxMinFlag);
DllExport int find_tangent_point(Point center, double R, Point Ps, 
					   Point Pm, Point V, Point* tangentPoint);
DllExport int expendShrinkSphere(Point* points, int n, Point* p, Sphere* maxMinSphere, int maxMinFlag);

DllExport Point* sphere_touchProbe_correction(Sphere sphere, Point* points, int n, 
											  double probeRadius, int getIn);

DllExport int sphere_bestfit_touchProbe( Point *points, int n, double tolerance,
         int iterateLimit, Sphere *init_sphere, int useInit, double *extraInfo,
         double probeRadius, int getIn, Point** correctedPoints);

DllExport int insidePyramid(Point* p, Point aPoint);

/*---------------------------------------------------------------*
 | CYLINDER --- implementation is in Cylinder.c                  |
 *---------------------------------------------------------------*/
DllExport double cylinder_sum_fun(double* point, double* b );
DllExport int cylinder_bestfit(Point* points, int n, double tolerance, int iterateLimit, 
                      Cylinder* init_cylinder, int useInit, double* extraInfo);
DllExport void cylinder_max_min(Cylinder cylinder, Point* points, int n, 
								Cylinder* maxMinCylinder, int maximum);
DllExport Point* cylinder_touchProbe_correction(Cylinder cylinder, Point* points, int n, 
												double probeRadius, int getIn);
DllExport int cylinder_bestfit_touchProbe( Point *points, int n, double tolerance,
         int iterateLimit, Cylinder   *init_cylinder, int    useInit,
         double *extraInfo, double probeR, int getIn, Point** correctedPoints);

/*---------------------------------------------------------------*
 | CONE --- implementation is in cone.c                          |
 *---------------------------------------------------------------*/
DllExport double cone_sum_fun(double* point, double* b );

DllExport int cone_bestfit(Point* points, int n, double tolerance, int iterateLimit, 
                        Cone* init_cone, int useInit, double* extraInfo);

DllExport Point* cone_touchProbe_correction(Cone cone, Point* points, int n, 
											 double probeRadius, int getIn);

DllExport int cone_bestfit_touchProbe( Point *points, int    n, double tolerance,
         int iterateLimit, Cone   *init_cone, int  useInit, double *extraInfo, 
         double probeRadius, int getIn, Point** correctedPoints);

DllExport int cone_bestfit2(Point* points, int n, double tolerance, int iterateLimit, 
                        Cone* init_cone, int useInit, double* extraInfo);

DllExport Point* cone_touchProbe_correction2(Cone cone, Point* points, int n, 
											 double probeRadius, int getIn);

DllExport int cone_bestfit_touchProbe2( Point *points, int    n, double tolerance,
         int iterateLimit, Cone *init_cone, int    useInit, double *extraInfo, 
         double probeRadius, int getIn, Point** correctedPoints);

/*---------------------------------------------------------------*
 | ERROR --- implementation is in error.c                        |
 *---------------------------------------------------------------*/
DllExport int error_coding(int error_code );
DllExport int error_check(char* location);

/*---------------------------------------------------------------*
 | SAMPLES --- implementation is in samples.c                    |
 *---------------------------------------------------------------*/
DllExport int  time_count(/*type, second, u_second */) ;
DllExport Point* generate_circle_points(Circle circle, int n, char* file_name );
DllExport Point* generate_ordered_circle_points(Circle circle, int n, char* file_name );
DllExport Point* generate_plane_points(Plane plane, int n, char* file_name, double dFact);
DllExport Point* generate_line_points(Line line, int n, char* file_name, double dFact);
DllExport Point* generate_sphere_points(Sphere sphere, int n, char* file_name );
DllExport Point* generate_cone_points(Cone cone, int n, char* file_name, double dFact );
DllExport Point* generate_ordered_cone_points(Cone cone, int n, char* file_name, double dFact );
DllExport Point* generate_cylinder_points(Cylinder cylinder, int n, char* file_name, double dFact );
DllExport Point* generate_ordered_cylinder_points(Cylinder cylinder, int n, char* file_name, double dFact);

DllExport Point* load_points(char* file_name, int* n);
DllExport void release_points(Point* points);

/*---------------------------------------------------------------*
 | STRUCTURE --- implementation is in structure.c                |
 *---------------------------------------------------------------*/
DllExport int  check_circle( Circle circle );
DllExport int  vector_2_circle(double* vector, Circle* circle );
DllExport int  circle_2_vector(Circle circle, double* vector );
DllExport void output_circle(Circle circle, char* message );

DllExport int  check_ellipse( Ellipse ellipse);
DllExport int  vector_2_ellipse(double* vector, Ellipse* ellipse);
DllExport int  ellipse_2_vector(Ellipse ellipse, double* vector );
//DllExport void output_ellipse(Ellipse ellipse, char* message );

DllExport int  check_arc( Arc arc );
DllExport int  vector_2_arc(double* vector, Arc* arc );
DllExport int  arc_2_vector(Arc arc, double* vector );
DllExport void output_arc(Arc arc, char* message );
 
DllExport int  check_plane( Plane plane );
DllExport int  vector_2_plane(double* vector, Plane* plane );
DllExport int  plane_2_vector( Plane plane, double* vector );
DllExport void output_plane(Plane plane, char* message );

DllExport int  check_line( Line line );
DllExport int  vector_2_line(double* vector, Line* line );
DllExport int  line_2_vector(Line line, double* vector );
DllExport void output_line(Line line, char* message );

DllExport int  check_sphere( Sphere sphere );
DllExport int  vector_2_sphere(double*  vector, Sphere* sphere );
DllExport int  sphere_2_vector( Sphere sphere, double* vector );
DllExport void output_sphere(Sphere sphere, char* message );

DllExport int  check_cone( Cone cone );
DllExport int  vector_2_cone( double* vector, Cone* cone );
DllExport int  cone_2_vector( Cone cone, double* vector );
DllExport void output_cone( Cone cone, char* message );

DllExport int  vector_2_cone2( double* vector, Cone* cone );
DllExport int  cone2_2_vector( Cone cone, double* vector );
 
DllExport int  check_cylinder(  Cylinder cylinder );
DllExport int  vector_2_cylinder( double* vector, Cylinder* cylinder );
DllExport int  cylinder_2_vector( Cylinder cylinder, double* vector );
DllExport void output_cylinder( Cylinder cylinder, char* message );
 
DllExport void output_vector(double* a, int n);
DllExport void output_triangle(double* a, int n );
DllExport void output_matrix(double* a, int n ) ;
 
/*---------------------------------------------------------------*
 | MARQUARDT --- implementation is in marquardt.c                |
 *---------------------------------------------------------------*/
DllExport int choleski_decompos(double* matrix, int row );
DllExport double* choleski_back(double* matrix, double* v, int row );
DllExport double* symmatrix_2_lowtriangle(double* matrix, double* low_triangle, int n );
DllExport double* lowtriangle_2_symmatrix(double* low_triangle, double* matrix, int n );
DllExport int marquardt(double* points, int m, double (*fun)(double*, double*), 
        double (*der_fun)(double*, double*, double*), 
        double* b, int n, double tolerance, int iterateLimit, double* sum_sq);

/*---------------------------------------------------------------*
 | 2D  GEOMETRY--- implementation is in geom_2d.c                |
 *---------------------------------------------------------------*/
DllExport double distsq_point_point_2D(Point point1, Point point2);
DllExport double distance_point_point_2D(Point point1, Point point2);

DllExport double distsq_circle_point_2D(Circle circle, Point point);
DllExport double distance_circle_point_2D(Circle, Point point);

DllExport double distsq_line_point_2D(Line line, Point point);
DllExport double distance_line_point_2D(Line line, Point point);

DllExport Point point_direction_move_2D(Point point, Point direction, double distance);

DllExport int intersection_line_line_2D(Line line1, Line line2, 
                              Point* intersection, double* angle);

DllExport int intersection_circle_line_2D(Circle circle, Line line, 
                                Point* p_intersect1, Point* p_intersect2);

DllExport int intersection_circle_circle_2D(Circle circle1, Circle circle2, 
                      Point* p_intersect1, Point* p_intersect2 );
DllExport int new_intersection_circle_circle(Circle circle1, Circle circle2,
							       Point *p_intersect1, Point *p_intersect2);

DllExport int intersection_arc_circle_2D(Arc arc, Circle circle, 
                      Point* p_intersect1, Point* p_intersect2 );

DllExport int intersection_arc_arc_2D(Arc arc1, Arc arc2, 
                      Point* p_intersect1, Point* p_intersect2 );

/*---------------------------------------------------------------*
 | TOLERANCE --- implementation is in toleranc.c                 |
 *---------------------------------------------------------------*/
DllExport double line_straightness_2D(Line line, Point* points, int n);
DllExport int line_profile_2D(Line line, Point* points, int n, 
                    double* profile_plus, double* profile_minus);

DllExport double circle_roundness_2D(Circle circle, Point* points, int n);
DllExport int circle_profile_2D(Circle circle, Point* points, int n, 
                    double* profile_plus, double* profile_minus);

DllExport double plane_flatness(Plane plane, Point* points, int n);
DllExport int plane_profile(Plane plane, Point* points, int n, 
                    double* profile_plus, double* profile_minus);

DllExport double sphere_roundness(Sphere sphere, Point* points, int n);
DllExport int sphere_profile(Sphere sphere, Point* points, int n,
                    double* profile_plus, double* profile_minus);

DllExport double line_per_ang_par( Point* points, int n, Point direction, int* index);

DllExport int cylinder_profile( Cylinder cylinder, Point* points, int n, 
                      double* profile_plus, double* profile_minus);

DllExport double cylinder_cylindricity( Cylinder cylinder, Point* points, int  n);

DllExport double truncate_decimals(double value);

DllExport void project2Plane(Point *coord, Point center, Point* points, Point* newPoint, int n);

/*---------------------------------------------------------------*
 | Width --- implementation is in Width.c					     |
 *---------------------------------------------------------------*/
DllExport int width_circle_point(Circle circle, Point point, WidthStruc* resWidth);
DllExport int width_line_point(Line line, Point point, WidthStruc* resWidth);
DllExport int width_point_point(Point point1, Point point2, WidthStruc* resWidth);
DllExport int width_line_line(Line line1, Line line2, WidthStruc* resWidth);
DllExport int width_line_circle(Line line, Circle circle, WidthStruc* resWidth);
DllExport int width_circle_circle(Circle circle1, Circle circle2, WidthStruc* resWidth);
DllExport int min_width_circle_point(Circle circle, Point point, WidthStruc* resWidth);
DllExport int min_width_line_point(Line line, Point point, WidthStruc* resWidth);
DllExport int min_width_point_point(Point point1, Point point2, WidthStruc* resWidth);
DllExport int min_width_line_line(Line line1, Point* line1_points, int num_line1, 
					Line line2, Point* line2_points, int num_line2,
					WidthStruc* resWidth);
DllExport int min_width_line_circle(Line line, Circle circle, WidthStruc* resWidth);
DllExport int min_width_circle_circle(Circle circle1, Circle circle2, WidthStruc* resWidth);
DllExport int max_width_circle_point(Circle circle, Point point, WidthStruc* resWidth);
DllExport int max_width_line_point(Line line, Point point, WidthStruc* resWidth);
DllExport int max_width_point_point(Point point1, Point point2, WidthStruc* resWidth);
DllExport int max_width_line_line(Line line1, Point* line1_points, int num_line1, 
					Line line2, Point* line2_points, int num_line2,
					WidthStruc* resWidth);
DllExport int max_width_line_circle(Line line, Circle circle, WidthStruc* resWidth);
DllExport int max_width_circle_circle(Circle circle1, Circle circle2, WidthStruc* resWidth);
DllExport int fcPerpendicularPoint(Line line, Point point, Point* perp_pt);
DllExport double fcCenterAngle(Line line1, Line line2);
DllExport double fcPerpendDist(double x, double y, double slope, double y_intercept); 
DllExport int GetAngle(Line line, double* return_angle);
DllExport double AngleBetweenInDegrees(Point origin, Point pt);
DllExport double AngleBetweenInRadians(Point origin, Point pt) ;

/*-------------------------------------------------------------------------*
 |  optimization fuction   ---- implemented in optim.c
 *-------------------------------------------------------------------------*/
DllExport int newton_minimize(double (*fun_der1)(double*, double), double (*fun_der2)(double*, double),
							  double *parameters, double *minimizer, double init_value, double tolerence);

DllExport double line_to_circle_sq_D1(double* param, double t);
DllExport double line_to_circle_sq_D2(double* param, double t);

DllExport double circle_to_circle_sq_D1(double* param, double t);
DllExport double circle_to_circle_sq_D2(double* param, double t);

/*---------------------------------------------------------------*
 | Contour --- implementation in Contour.c						 |
 *---------------------------------------------------------------*/
DllExport int contour_center_fit(Point* points, int n, Contour* init_contour);
DllExport int contour_min_x(Point* points, int n, Contour* init_contour);
DllExport int contour_max_x(Point* points, int n, Contour* init_contour);
DllExport int contour_min_y(Point* points, int n, Contour* init_contour);
DllExport int contour_max_y(Point* points, int n, Contour* init_contour);
DllExport int contour_min_z(Point* points, int n, Contour* init_contour);
DllExport int contour_max_z(Point* points, int n, Contour* init_contour);
DllExport int contour_min_angle(Point* points, int n, Contour* init_contour, double angle);
DllExport int contour_max_angle(Point* points, int n, Contour* init_contour, double angle);
DllExport int contour_max_radius(Point* points, int n, Contour* init_contour);
DllExport int contour_min_radius(Point* points, int n, Contour* init_contour);

#ifdef __cplusplus
}
#endif

#endif

