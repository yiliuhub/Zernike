
/******************************************************************************
 * FILE: distance.c
 *       Implementation of functions for distances between objects.
 *
 *       double distsq_circle_point()
 *       double distsq_line_line()
 *       double distsq_line_plane()
 *       double distsq_line_point()
 *       double distsq_plane_plane()
 *       double distsq_plane_point()
 *       double distsq_point_point()
 *       double distsq_point_sphere()
 *       double distsq_cone_point()
 *       double distsq_cylinder_point()
 *       
 *       double distance_circle_point()
 *       double distance_line_line()
 *       double distance_line_plane()
 *       double distance_line_point()
 *       double distance_plane_plane()
 *       double distance_plane_point()
 *       double distance_point_point()
 *       double distance_point_sphere()
 *       double distance_cone_point()
 *       double distance_cylinder_point()
 *        
 * 10/05/95  LXD  Initial creation.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "geometry.h"

/******************************************************************************
 *     POINT and POINT
 *     point1 (X1, Y1, Z1) 
 *     point2 (X2, Y2, Z2) 
 *     distance square   d^2 = (X1-X2)^2 + (Y1-Y2)^2 + (Z1-Z1)^2 
 ******************************************************************************/
double
distsq_point_point( Point point1, Point point2)
{
   double dist_x, dist_y, dist_z;/* distances in each direction to middle point */
   double distance_square;       /* square of the distance */
     
     
   dist_x = point1.x - point2.x;
   dist_y = point1.y - point2.y;
   dist_z = point1.z - point2.z;
     
   distance_square = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
 
   return distance_square;
}

double
distance_point_point( Point point1, Point point2)
{
   return sqrt(distsq_point_point(point1, point2));
}

/****************************************************************************** 
 *     POINT and LINE 
 *     point (X1, Y1, Z1)
 *     Line through point (X2,Y2,Z2) and direction is (L,M,N) 
 *     square of distance
 *     d^2 = (X1-X2)^2 + (Y1-Y2)^2 + (Z1-Z1)^2  
 *           - [L*(X1-X2) + M*(Y1-Y2) + N*(Z1-Z2)]^2/(L^2 + M^2 + N^2)
 *
 ******************************************************************************/
 
double
distsq_line_point( Line  line, Point point)
{
   double dist_x, dist_y, dist_z;/* distances in each direction to middle point */
   double lx_my_nz;              /* temp variable to save the value for square */
   double l, m, n;               /* direction of the line */
   double distance_square;       /* square of the distance */
   

   dist_x = point.x - line.middle_point.x;
   dist_y = point.y - line.middle_point.y;
   dist_z = point.z - line.middle_point.z;

   l = line.direction.x;
   m = line.direction.y;
   n = line.direction.z;

   lx_my_nz = l*dist_x + m*dist_y + n*dist_z;  

	if (0.0 != l*l+m*m+n*n)
	{
		distance_square = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z
              - lx_my_nz*lx_my_nz/(l*l+m*m+n*n);
	}
	else
	{
		distance_square = 0.0;
	}

   if (distance_square <= 0.0)
      distance_square = 0.0;

   return distance_square;
}

double
distance_line_point( Line  line,   Point point)
{

   return sqrt(distsq_line_point(line, point));
}
/******************************************************************************  
 *    LINE and LINE 
 *    Line 1 : through (X1,Y1,Z1) and direction is (L1, M1, N1)
 *    Line 2 : through (X2,Y2,Z2) and direction is (L2, M2, N2)
 *       case 1 (L1, M1, N1) = t(L2, M2, N2)
 *              Two lines are parallel.
 *              d = distance from point (X1,Y1,Z1) to Line 2 
 *
 *       case 2 (L1, M1, N1) != (L2, M2, N2)
 *              distance:  d = |d'| where
 *
 *                     |M1 N1|            |N1 L1|            |L1 M1|
 *            (X1 - X2)|M2 N2| + (Y1 - Y2)|N2 L2| + (Z1 - Z2)|L2 M2|
 *     d' =  ________________________________________________________ 
 *                   ___________________________________
 *                  /        2           2           2
 *                 /  |M1 N1|  +  |N1 L1|  +  |L1 M1|
 *               \/   |M2 N2|     |N2 L2|  +  |L2 M2|
 ******************************************************************************/ 
  
double 
distance_line_line( Line  line1, Line  line2)
{ 
   double l1, m1, n1; /* stand for L1, M1 and N1 in above formula        */
   double l2, m2, n2; /* stand for L2, M2 and N2 in above formula        */
   double det_mn;     /* result of determinants to save repeat calulates */
   double det_nl;     /*        |M1 N1|        |N1 L1|          |L1 M1|  */
   double det_lm;     /* det_mn=|M2 N2| det_nl=|N2 L2| det_lm = |L2 M2|  */
   double numerator;  /* result of above fornula's numerator             */
   double denominator;/* result of above fornula's denominator           */
   double distance;   /* distance = numerator/denominator                */

   l1 = line1.direction.x ;
   m1 = line1.direction.y ;
   n1 = line1.direction.z ;

   l2 = line2.direction.x ;
   m2 = line2.direction.y ;
   n2 = line2.direction.z ;

   det_mn = m1*n2 - n1*m2 ; 
   det_nl = n1*l2 - l1*n2 ; 
   det_lm = l1*m2 - m1*l2 ; 

   if (det_mn==0 && det_nl==0 && det_lm==0) 
      /* case 1: two lines are parallel */
      distance = distance_line_point(line2, line1.middle_point); 
   else /* case 2: two lines are not parallel */
   {
      numerator = (line1.middle_point.x - line2.middle_point.x )*det_mn
                + (line1.middle_point.y - line2.middle_point.y )*det_nl
                + (line1.middle_point.z - line2.middle_point.z )*det_lm ;
      denominator = sqrt(det_mn*det_mn + det_nl*det_nl + det_lm*det_lm);
	  if (0.0 != denominator)
	  {
		  distance = fabs(numerator/denominator);
	  }
	  else
	  {
		  distance = 0.0;
	  }
   }
   return distance;
} 

double 
distsq_line_line( Line  line1, Line  line2)
{
   double distance;

   distance = distance_line_line(line1, line2);

   return distance*distance;
}
/******************************************************************************
 *     POINT and PLANE 
 *     Point (X0, Y0, Z0)
 *     Plane Through (X1, Y1, Z1) and direction of normal is (A, B, C)
 *
 *     the distance from the Point to the Plane
 *
 *               |A1(X0-X1) + B1(Y0-Y1) + C1(Z0-Z1)|
 *        d =    ___________________________________
 *                       ________________
 *                      /   2    2    2
 *                    \/  A1 + B1 + C1
 ******************************************************************************/
double
distsq_plane_point( Plane plane, Point point)
{
    double dist = distance_plane_point(plane, point);

    return dist*dist;
}

double
distance_plane_point( Plane plane, Point point)
{
    double a1, b1, c1;  /* stand for A1, B1 and C1 in above formula */
    double numerator;   /* result of numerator in above formula     */
    double denominator; /* result of denominator in above formula   */ 
    double dist;    /* distance = numerator/denominator         */

    a1 = plane.direction.x;
    b1 = plane.direction.y;
    c1 = plane.direction.z;

    numerator = a1 * (point.x - plane.center.x)
              + b1 * (point.y - plane.center.y)
              + c1 * (point.z - plane.center.z);
    denominator =  sqrt(a1*a1 + b1*b1 + c1*c1) ;
	if (0.0 != denominator)
	{
		dist = numerator/denominator;
	}
	else
	{
		dist = 0.0;
	}

    return dist;
}

/******************************************************************************
 *    LINE and PLANE
 *    Line through Point (X1, Y1, Z1) in direction (L1, M1, N1)
 *    Plane through Point (X2, Y2, Z2) in (normal) direction (A2, B2, C2)
 *
 *    case 1: L1*A2 + M1*B2 + N1*C2 != 0
 *            The Line is not parallel to the Plane distance d = 0.
 *    case 2: L1*A2 + M1*B2 + N1*C2 == 0
 *            The Line is parallel to the Plane
 *            d = distance from the line's middle_point to the plane.
 ******************************************************************************/

double 
distance_line_plane( Line  line, Plane plane)
{
   double l1, m1, n1; /* stand for L1, M1, N1 in above formula */
   double a2, b2, c2; /* stand for A2, B2, C2 in above formula */
   double distance;   /* distance from the line to the plane   */

   l1 = line.direction.x;
   m1 = line.direction.y;
   n1 = line.direction.z;

   a2 = plane.direction.x;
   b2 = plane.direction.y;
   c2 = plane.direction.z;

   if ((l1*a2 + m1*b2 + n1*c2) != 0) 
      /* case 1: The Line and the Plane are intersected */
      distance = 0.0;
   else /* case 2 The Line and the Plane are parallel */
      distance = distance_plane_point(plane, line.middle_point);

   return distance;
}
 
double
distsq_line_plane( Line  line, Plane plane)
{
   double distance;
   
   distance = distance_line_plane(line, plane);
   return distance*distance;
}

/****************************************************************************** 
 *    PLANE and PLANE 
 *   
 *    Plane 1 through (X1, Y1, Z1) in direction	(A1, B1, C1)
 *    Plane 2 through (X2, Y2, Z2) in direction	(A2, B2, C2)
 *   
 *     case 1: Plane 1 || Plane2
 *             distance is from (X1, Y1, Z1) to Plane 2)
 *     case 2: Plane 1 intersects Plane 2 
 *             distance is 0
 ******************************************************************************/ 

double  
distance_plane_plane( Plane plane1, Plane plane2 )
{ 
   double distance = 0.0;

   if (parallel_plane_plane(plane1, plane2))
      distance = distance_plane_point(plane2, plane1.center);
   
   return distance;
} 
  
/******************************************************************************
 *    POINT and CIRCLE
 *
 *    Point  P=(X1, Y1, Z1)
 *    Circle 2 center O=(X2, Y2, Z2) in direction (A2, B2, C2) with radius r
 *        d 
 *    P.______|
 *    /|h     |
 *  Q/_|______|_______
 *     P' r   O   r
 *
 *       2      2       2       2   2         2
 *   dist = |PQ| = |PP'| + |P'Q| = h + (r - d)
 *
 ******************************************************************************/
 
double
distance_circle_point( Circle circle, Point  point)
{ 
   double vecCircle[7];

   double distance ;
 
   circle_2_vector(circle, vecCircle);
   distance = circle_sum_fun((double*)(&point), vecCircle);

   return distance;
} 

double
distsq_circle_point( Circle circle, Point  point)
{ 
   double distance ;
 
   distance = distance_circle_point(circle, point); 
 
   return distance*distance;
}


/******************************************************************************
 *    POINT and SPHERE
 *
 *    Point  P=(X1, Y1, Z1)
 *    Sphere 2 center O=(X2, Y2, Z2)  radius r
 *
 *   dist = |PO| - r
 *
 ******************************************************************************/
 
double
distance_point_sphere( Point  point, Sphere sphere)
{
   double vecSphere[4];

   double distance ;

   sphere_2_vector(sphere, vecSphere);
   distance = sphere_sum_fun((double*)(&point), vecSphere);

   return distance;
}

double
distsq_point_sphere( Point point, Sphere sphere)
{
   double distance ;
 
   distance = distance_point_sphere(point, sphere);
 
   return distance*distance;
}
 

/******************************************************************************
 *    POINT and CONE
 *
 ******************************************************************************/
 
double
distance_cone_point( Cone cone, Point  point)
{
	double distance;
	double vecCone[7];
 
	cone_2_vector(cone, vecCone);
	distance = cone_sum_fun((double*)(&point), vecCone);
 
/*
	Point VToP;
	Line cone_axis;
	double height;
	double difference;
	double projected_length;
	////////////////////////////////////////////////////////////
	//	The following code calculates a non-negative distance 
	//	from a point to a cone.
	// Yi Liu 11/99
	////////////////////////////////////////////////////////////


	//vertex to point
	VToP = v_sub(point, cone.vertex);
				
	cone_axis.middle_point = cone.vertex;
	cone_axis.direction = cone.direction;
	cone_axis.direction = unit_direction(cone_axis.direction);

	distance = distance_line_point(cone_axis, point);
	projected_length = fabs(v_dot(VToP, cone_axis.direction));
	height = projected_length * tan(cone.angle/2);

	distance = fabs(distance - height) * cos(cone.angle/2);
*/
   return distance;
}
 
double
distsq_cone_point(Cone cone, Point point)
{
   double distance ;
 
   distance = distance_cone_point(cone, point);
 
   return distance*distance;
}
 
/******************************************************************************
 *    POINT and CYLINDER
 *
 ******************************************************************************/
 
double
distance_cylinder_point( Cylinder cylinder, Point  point)
{
   double vecCylinder[7];
 
   double distance ;
 
   cylinder_2_vector(cylinder, vecCylinder);
   distance = cylinder_sum_fun((double*)(&point), vecCylinder);
 
   return distance;
}
 
double
distsq_cylinder_point(Cylinder cylinder, Point point)
{
   double distance ;
 
   distance = distance_cylinder_point(cylinder, point);
 
   return distance*distance;
}
 


