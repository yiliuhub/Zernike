#include "stdafx.h"
#include <math.h>
#include <stdlib.h>

#include "marquard.h"



Point v_create(double x, double y, double z)
{
	Point p;
	p.x = x, p.y = y, p.z = z;
	return p;
}

double distance_point_point(Point p, Point q)
{
	return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y) + (p.z-q.z)*(p.z-q.z));
}
