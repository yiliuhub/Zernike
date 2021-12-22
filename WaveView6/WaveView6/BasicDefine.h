#pragma once

namespace LLW
{
	template< typename T >
	struct Point2D
	{
		T x;
		T y;
	};

	template< typename T >
	struct Point3D
	{
		T x;
		T y;
		T Z;
	};

	typedef		Point2D < char >			Point2c;

	typedef		Point2D < unsigned char >	Point2u;

	typedef		Point2D < int >				Point2n;

	typedef		Point2D < float >			Point2f;
	
	typedef		Point2D < double >			Point2d;

	typedef		Point3D < char >			Point3c;

	typedef		Point3D < unsigned char >	Point3u;
	
	typedef		Point3D < int >				Point3n;
	
	typedef		Point3D < float >			Point3f;
	
	typedef		Point3D < double >			Point3d;
}