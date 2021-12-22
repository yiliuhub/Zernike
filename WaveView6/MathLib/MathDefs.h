
#ifndef MATH_BASIC_H
#define MATH_BASIC_H

#ifndef PI
#define PI					3.1415926535897932384626433832795
//							3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647
#endif

#ifndef PIBYTWO
  #define PIBYTWO				1.5707963267948966192313216916398
#endif

#define SMALL_REAL			0.0000000001
#define SMALL_ANGLE			0.00001
#define DIVISOR_CHECK		0.000001
#define	EPS					1.e-100
#define INCH_MM_CONVERSION	25.4


typedef struct _Xyz_Vector {
	double x;		/* x y z coordinate set */
	double y;
	double z;
} Xyz_Vector;

typedef struct _Matrix {
	Xyz_Vector a;
	Xyz_Vector b;
	Xyz_Vector c;	
} Matrix;

typedef struct _Transform {
	Xyz_Vector d;	/* 3x1 translation */
	Matrix m;		/* 3x3 transform */
} Transform;

#endif