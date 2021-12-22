/**************************************************************
*
*	llwZernike.cpp 
*	Author:	Yi Liu
*	Date:	May 2007
*
****************************************************************/
#include "StdAfx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "..\MathLib\vmblock.h"
#include "..\MathLib\u_proto.h"
#include "LLWzernike.h"

#include "LLW_basics.h"


#define EPS 0.000000000000000000000000001


LLWzernike::LLWzernike():
m_dF(1),
m_nSize(0),
m_nOrder(ZERNIKE_ORDER),
m_pRefX(NULL),
m_pRefY(NULL),
m_pDeltaX(NULL),
m_pDeltaY(NULL),
//m_error(NO_ERROR),
	m_ppMatrix(NULL),
	m_pB(NULL),
	m_pX(NULL),
	m_dScalor(1)
{
}

LLWzernike::~LLWzernike()
{
	Reset();
}
void LLWzernike::Reset()
{
	//m_error = NO_ERROR;
	m_nSize = 0;
	m_dF = 1.0;
	
	delete [] m_pRefX;
	delete [] m_pRefY;
	delete [] m_pDeltaX;
	delete [] m_pDeltaY;
	m_pRefX = NULL;
	m_pRefY = NULL;
	m_pDeltaX = NULL;
	m_pDeltaY = NULL;
}

bool LLWzernike::ValidityCheck()
{
	bool bOk = true;

	if(	(m_pRefX == NULL) || (m_pRefY == NULL) ||
		(m_pDeltaX == NULL)|| (m_pDeltaY == NULL) )
	{
		bOk = false;
	}

	return bOk;
}

/******************************************************************************************
 SetPointList(unsigned short size, double* pX, double* pY, double* pDraftX, double* pDraftY)
 Sets Reference point (pX, pY) and its measured points (pDraftX, pDraftY).
 The delta will be ( pDraftX - pX, pDraftY - pY )
******************************************************************************************/
bool LLWzernike::SetPointList(unsigned short size, double* pX, double* pY, double* pDraftX, double* pDraftY)
{
	bool bOk = true;
	Reset();

	// insufficient data
	if(	(2*size < 0.5*(ZERNIKE_ORDER+1)*(ZERNIKE_ORDER+2)) ||
		(pX == NULL) || (pY == NULL) || (pDraftX == NULL) || (pDraftY == NULL))
	{
		bOk = false;
	}

	if(bOk)
	{
		m_nSize = size;
		m_pRefX = new double[size];
		m_pRefY = new double[size];
		m_pDeltaX = new double[size];
		m_pDeltaY = new double[size];

		if( (m_pRefX == NULL)	|| (m_pRefX == NULL) ||
			(m_pDeltaX == NULL) || (m_pDeltaY == NULL))
		{
			//m_error = MEMORY_ALLOCATION_FAILURE;
			bOk = false;
		}
	}

	if(bOk)
	{
		for(int i=0; i<size; i++)
		{
			m_pRefX[i] = pX[i];
			m_pRefY[i] = pX[i];
			m_pDeltaX[i] = pDraftX[i] - pX[i];
			m_pDeltaY[i] = pDraftY[i] - pY[i];
		}
	}

	if(!bOk)
	{
		Reset();
	}

	return bOk;
}

bool LLWzernike::LoadDataFromFile()
{
	int index = 0;
	int count = 0;
	char str[100];
	double x;
	double f[5];
	bool bOk = true;

	FILE* pData = fopen("C:\\My Document\\data.txt", "r");

	if(pData == NULL)
	{
		bOk = false;
		//m_error = FAILED_OPEN_FILE;
	}

	Reset();
	m_nSize = 1000;
	m_pRefX = new double[m_nSize];
	m_pRefY = new double[m_nSize];
	m_pDeltaX = new double[m_nSize];
	m_pDeltaY = new double[m_nSize];

	if(m_pRefX == NULL || m_pRefY == NULL || m_pDeltaX == NULL || m_pDeltaY == NULL)
	{
		bOk = false;
		//m_error = MEMORY_ALLOCATION_FAILURE;
		Reset();
	}

	if(bOk)
	{
		m_nSize = 0;
		rewind(pData);
		while(!feof(pData))
		{
			count = fscanf (pData, "%s", str);
			x = atof(str);
			
			if(count > 0)
			{
				if(x==1.0)
				{
					m_pRefX[m_nSize] = f[(index-4+4)%4];
					m_pRefY[m_nSize] = f[(index-3+4)%4];
					m_pDeltaX[m_nSize] = (f[(index-2+4)%4] - f[(index-4+4)%4]);
					m_pDeltaY[m_nSize] = (f[(index-1+4)%4] - f[(index-3+4)%4]);
					m_nSize++;
				}
				else
				{
					f[index] = x;
					index = (index+1)%4;
				}
			}
		}
	}

	fclose(pData);

	FILE* pReadData = fopen("C:\\My Document\\read_data.txt", "w");

	for(int i=0; i<m_nSize; i++)
	{
		fprintf(pReadData, "%20.12f  %20.12f  %20.12f  %20.12f\n", m_pRefX[i], m_pRefY[i],  m_pRefX[i]+m_pDeltaX[i],  m_pRefY[i]+m_pDeltaY[i]);
	}
	
	fclose(pReadData);

	return bOk;
}

bool LLWzernike::PreProcessData()
{
	bool bOk = true;
	int i = 0;
	double ctrX, ctrY;
	double minX, maxX, minY, maxY;

	bOk = ValidityCheck();
	if(bOk)
	{
		minX = maxX = m_pRefX[0];
		minY = maxY = m_pRefY[0];

		// seeking min/max
		for(i= 0; i<m_nSize; i++)
		{
			minX = (minX < m_pRefX[i])? minX : m_pRefX[i];
			minY = (minY < m_pRefY[i])? minY : m_pRefY[i];
			maxX = (maxX > m_pRefX[i])? maxX : m_pRefX[i];
			maxY = (maxY > m_pRefY[i])? maxY : m_pRefY[i];
		}
		// set center
		ctrX = 0.5*(minX + maxX);
		ctrY = 0.5*(minY + maxY);
	}
	// normalizing data, we use focus length to normalize data currently
	m_dScalor = 0.5*__max( maxX-minX, maxY-minY );
	m_dScalor = __max(m_dScalor, 1.0);

	for(i=0; i<m_nSize; i++)
	{
		m_pRefX[i] = (m_pRefX[i] - ctrX)/m_dScalor;
		m_pRefY[i] = (m_pRefY[i] - ctrY)/m_dScalor;
	//	m_pDeltaX[i] /= m_dF;
	//	m_pDeltaY[i] /= m_dF;
	}

	return bOk;
}

bool LLWzernike::MapIJtoK(int i, int j, int& k)
{
	bool bOk = true;
	if(i > ZERNIKE_ORDER || i< 0 || i < j)
	{
		bOk = false;
	}
	if(bOk)
	{
		k = (int)(0.5*(i)*(i+1) + j);
	}

	return bOk;
}

/////////////////////////////////////////////////////////////////////
// k = 0.5*i*(i+1) + j  with j <= i
// so 2K >= i(i+1)
//  (i+0.5)^2 <= 2K + 0.25
//	i <= sqrt(2k + 0.25) - 0.5
/////////////////////////////////////////////////////////////////////
bool LLWzernike::MapKtoIJ(int k, int& i, int& j)
{
	bool bOk = true;
	if((k < 0) || (k >= 0.5*(ZERNIKE_ORDER+1)*(ZERNIKE_ORDER+2)))
	{
		bOk = false;
	}
	if(bOk)
	{
		i = (int)(sqrt(2.*k+0.25) - 0.5);
		j = (int)(k - 0.5*(i)*(i+1));
	}
	return bOk;
}

/*******************************************************************************
	let A*w = b
	then the solution for the minimizing of the linear problem is:
	 w = (A'A)^(-1)A'b
	 here A' represent the transpose matrix, ^(-1) represents the inverse matrix
********************************************************************************/
bool LLWzernike::SetLinearSystem()
{
	bool bOk = true;
	int i,j;
	// total number of zernike polynomials
	m_nOrder = 8;
	int count = (int)(0.5*((m_nOrder+1)*(m_nOrder+2)));

	int row = 2*m_nSize+1, 
		col = count;


	double* b = new double[row];
	double* w = new double[col];
	memset(w, 0, col);

	// allocate row by col matrix, this is the coefficient matrix in Zernike System
	double** A = new double* [row];
	for (i = 0 ; i < row; i++)
		A[i] = new double[col];

	// B will be the transpose of A
	double** B = new double* [col];
		for (i = 0 ; i < col ; i++)
			B[i] = new double[row];

		// C is a square matrix
	double** C = new double* [col];
		for (i = 0 ; i < col ; i++)
			C[i] = new double[col];

	// Bb is the vector B*b
	double* Bb = new double[col];

	int m, n;
	//setting matrix A
	for(i=0; i< m_nSize; i++)
	{
		// this is the c term in the document by Patrick Y. Maeda
		b[i] = m_pDeltaX[i]/m_dF;
		b[i+m_nSize] = m_pDeltaY[i]/m_dF;

		for(j=0; j<col;j++)
		{
			MapKtoIJ(j, m, n);
			// this is the h term in the document by Patrick Y. Maeda
			A[i][j]			= ZernikePoly_Dx(m, n, m_pRefX[i], m_pRefY[i]);
			A[i+m_nSize][j] = ZernikePoly_Dy(m, n, m_pRefX[i], m_pRefY[i]);
		}
	}

	// last equation, 为什么系数和为零？
	b[row-1] = 100;
	for(j=0; j<col;j++)
	{	A[2*m_nSize][j] = 1;
	}

	// compute transpose matrix B of A。 B为A的转置矩阵
	matrix_tran(row, col, A, B);
	// compute B*A
	matrix_prod(col, row, col, B, A, C);
	// compute B*b
	matrix_prod_vec(col, row, B, b, Bb);

	// solve linear equation (A'A)w = A'b
	// here A' = B, the transpose of A
	int result = 0;
	int signd = 1;
	int* perm = new int[count];
	double** lumat = new double* [count];
		for (i = 0 ; i < count ; i++)
			lumat[i] = new double[count];

	REAL determinant = det(col,C);
	if(fabs(determinant) > 0.0000000000001)
		result = gauss(0, col, C, lumat, perm, Bb, w, &signd);
	else
		bOk = false;

	//dump data
	FILE* matrixC = fopen("C:\\My Document\\matrixC.txt", "w");

	// reference points and delta
	fprintf(matrixC, "\n\n RefY  RefX  DeltaY DeltaX \n");
	for(i=0; i<m_nSize;i++)
	{	
		fprintf(matrixC, "%12.8f %12.8f %12.8f %12.8f \n", m_pRefY[i], m_pRefX[i], m_pDeltaY[i], m_pDeltaX[i]);
	}

	// matrix of the linear system
	fprintf(matrixC, "\n\n matrix A=\n");
	for(i=0; i<row;i++)
	{	for(j=0;j<col;j++)
		{	fprintf(matrixC, "%12.5f ", A[i][j]);
		}
		fprintf(matrixC, "\n");
	}


	// matrix of the linear system
	fprintf(matrixC, "\n\n matrix C=\n");
	for(i=0; i<col;i++)
	{	
		for(j=0;j<col;j++)
		{	fprintf(matrixC, "%12.5f ", C[i][j]);
		}
		fprintf(matrixC, "\n");
	}

	// non-homogeneous term
	fprintf(matrixC, "\n b= \n");
	for(j=0;j<row;j++)
	{	fprintf(matrixC, "%12.8f\n ", b[j]);
	}

	// solution of Zernike Coefficients
	fprintf(matrixC, "\n w: \n");
	for(j=0;j<col;j++)
	{
		MapKtoIJ(j, m,n);
		w[j] = m_dScalor*PIXEL_MICRON_RATIO*w[j];
		fprintf(matrixC, "w(%d, %d) = %12.8f\n", m,n, w[j]);
	}
	fclose(matrixC);
	//end of dumping

	//release memory
	for(i=0; i<row;i++)
		delete[] A[i];
	for(i=0; i<col;i++)
		delete[] B[i];
	for(i=0; i<col;i++)
		delete[] C[i];
	for(i=0; i<count;i++)
		delete[] lumat[i];

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] lumat;

	delete[] Bb;
	delete[] b;
	delete[] w;
	delete[] perm;
	
	return bOk;
}

bool LLWzernike::GoZernike()
{
	bool bOk = true;
	TestDataIOFunction();

	bOk = LoadDataFromFile();
	if(bOk)
	{
		bOk = PreProcessData();
		if(bOk)
		{
			bOk = SetOrder(8);
			bOk = SetFocusLength(700);
			bOk = SetLinearSystem();
		}
	}

	return bOk;
}

bool LLWzernike::SetFocusLength(double focus)
{
	bool bOk = focus>0? true:false;
	if(bOk) 
		m_dF = focus;
	else
		//m_error = INVALID_FOCUS_LENGTH;

	return bOk;
}

bool LLWzernike::SetNormalScalor(double scalor)
{
	bool bOk = true;
	if(scalor <=0)
		bOk = false;
	else
		m_dScalor = scalor;
	return bOk;
}

bool LLWzernike::SetOrder(int order)
{
	bool bOk = order>0? true:false;
	if(bOk)
	{	if(order > ZERNIKE_ORDER)
			m_nOrder = ZERNIKE_ORDER;
		else
			m_nOrder = order;
	}
	else
	{	
		//m_error = INVALID_ORDER;
	}

	return bOk;
}

double LLWzernike::Zernike(double x, double y)
{
	int i, j;
	double z = 0.0;

	for(i=0; i<=ZERNIKE_ORDER; i++)
	{
		for(j=0; j<=i; j++)
		{
			z += m_dCoefficient[i][j]*ZernikePoly(i, j, x, y);
		}
	}
	return z;
}

void LLWzernike::TestDataIOFunction()
{
	FILE* polydata = fopen("C:\\My Document\\Zernike_Poly_Data.txt", "w");

	int size = (int)(0.5*(m_nOrder+1)*(m_nOrder+2));
	double x = 3.1, y = 2.7;

	int m, n;
	for(int i= 0; i<size; i++)
	{
		MapKtoIJ(i, m, n);
		fprintf(polydata, "%20.10f   %20.10f   %20.10f\n", ZernikePoly(m,n,x,y), ZernikePoly_Dx(m,n,x,y), ZernikePoly_Dy(m,n,x,y));
	}
	fclose(polydata);
  return;
}

double LLWzernike::Poly(double a, double b, double x, double y)
{
	return pow(x,a)*pow(y,b);
}
double LLWzernike::Poly_Dx(double a, double b, double x, double y)
{
	if(a<1 && fabs(x)<EPS)
		return 0.0;
	else
		return (a*Poly(a-1,b, x, y));
}
double LLWzernike::Poly_Dy(double a, double b, double x, double y)
{
	if(b<1 && fabs(y)<EPS)
		return 0.0;
	else
		return (b*Poly(a,b-1, x, y));
}

/////////////////////////////////////////////////////////////////
//  This is the Zernike Polynomila of Z(m, n)
//
/////////////////////////////////////////////////////////////////
double LLWzernike::ZernikePoly(unsigned int m, unsigned int n, double x, double y)
{
		if((m>ZERNIKE_ORDER) || (n>ZERNIKE_ORDER))
	{
		return 0.0;
	}

	double Z_mn = 0.0;
	int i=0, j=0;

	switch(m)
	{
	case 0:	 //Z(0,*)
			Z_mn = Poly(0,0,x, y);
		break;	
	case 1:	 //Z(1,*)
		switch(n){
		case 0:	Z_mn = Poly(0,1, x, y);
			break;	
		case 1: Z_mn = Poly(1,0, x, y);
			break;
		}
		break;	
	case 2:	 //Z(2,*)
		switch(n){
		case 0: Z_mn = (2*Poly(1,1,x,y));	//Z(2,0) = 2xy
				break;	
		case 1: Z_mn = (-Poly(0,0,x,y)+2*Poly(0,2,x,y)+2*Poly(2,0,x,y));	//Z(2,1) = 1+2x2 + 2y2
				break;	
		case 2: Z_mn = (-Poly(0,2,x,y)+Poly(2,0,x,y));	//Z(2,2) = -x2 + y2
			break;
		}
		break;
	case 3:	 //Z(3,*)
		switch(n){
		case 0: Z_mn = (-Poly(0,3,x,y)+3*Poly(2,1,x,y));// Z(3,0) = -x3 + 3xy2
				break;	
		case 1: Z_mn = (-2*Poly(0,1,x,y)+3*Poly(0,3,x,y)+3*Poly(2,1,x,y));	// Z(3,1) = -2x + 3x3 + 3xy2 ((0),(-2,0),(0,0,0),(3,0,3,0))
				break;	
		case 2: Z_mn = (-2*Poly(1,0,x,y)+3*Poly(3,0,x,y)+3*Poly(1,2,x,y));// Z(3,2) = -2y + 3y3 + 3x2y ((0),(0,-2),(0,0,0),(0,3,0,3))
				break;	
		case 3: Z_mn = (Poly(3,0,x,y)-3*Poly(1,2,x,y));// Z(3,3) = y3 - 3x2y ((0),(0,0),(0,0,0),(0,-3,0,1))
				break;
		}
		break;
	case 4:	 //Z(4,*)
		switch(n){
		case 0: Z_mn = (-4*Poly(1,3,x,y)+4*Poly(3,1,x,y));// Z(4,0) = -4x3y + 4xy3 ((0),(0,0),(0,0,0),(0,0,0,0),(0,-4,0,4,0))
				break;	
		case 1: Z_mn = (-6*Poly(1,1,x,y)+8*Poly(1,3,x,y)+8*Poly(3,1,x,y));// Z(4,1) = -6xy + 8x3y + 8xy3 ((0),(0,0),(0,-6,0),(0,0,0,0),(0,8,0,8,0))
				break;	
		case 2: Z_mn = (Poly(0,0,x,y)-6*Poly(0,2,x,y)-6*Poly(2,0,x,y)+6*Poly(0,4,x,y)+12*Poly(2,2,x,y)+6*Poly(4,0,x,y));// Z(4,2) = 1 - 6x2 - 6y2 + 6x4 + 12x2y2 + 6y4 ((1),(0,0),(-6,0,-6),(0,0,0,0),(6,0,12,0,6))
				break;	
		case 3: Z_mn = (3*Poly(0,2,x,y)-3*Poly(2,0,x,y)-4*Poly(0,4,x,y)+4*Poly(4,0,x,y));// Z(4,3) = 3x2 - 3y2 - 4x4 + 4y4 ((0),(0,0),(3,0,-3),(0,0,0,0),(-4,0,0,0,4))
				break;	
		case 4: Z_mn = (Poly(0,4,x,y)-6*Poly(2,2,x,y)+Poly(4,0,x,y));// Z(4,4) = x4 - 6x2y2 + y4 ((0),(0,0),(0,0,0),(0,0,0,0),(1,0,-6,0,1))
				break;
		}
		break;
	case 5:	 //Z(5,*)
		switch(n){
		case 0: Z_mn = (   Poly(0,5,x,y)-10*Poly(2,3,x,y)+5*Poly(4,1,x,y));// Z(5,0) = x5 - 10x3y2 + 5xy4
				break;	
		case 1: Z_mn = ( 4*Poly(0,3,x,y)-12*Poly(2,1,x,y)-5*Poly(0,5,x,y)+10*Poly(2,3,x,y)+15*Poly(4,1,x,y));//Z(5,1) = 4x3 - 12xy2 - 5x5 + 10x3y2 + 15xy4
				break;
		case 2: Z_mn = ( 3*Poly(0,1,x,y)-12*Poly(0,3,x,y)-12*Poly(2,1,x,y)+10*Poly(0,5,x,y) + 20*Poly(2,3,x,y) + 10*Poly(4,1,x,y));//Z(5,2) = 3x - 12x3 - 12xy2 + 10x5 + 20x3y2 + 10xy4
				break;	
		case 3: Z_mn = ( 3*Poly(1,0,x,y)-12*Poly(3,0,x,y)-12*Poly(1,2,x,y)+10*Poly(5,0,x,y)+20*Poly(3,2,x,y)+10*Poly(1,4,x,y));//Z(5,3) = 3y - 12y3 - 12x2y + 10y5 + 20x2y3 + 10x4y
				break;	
		case 4: Z_mn = (-4*Poly(3,0,x,y)+12*Poly(1,2,x,y)+5*Poly(5,0,x,y)-10*Poly(3,2,x,y)-15*Poly(1,4,x,y));//Z(5,4) = -4y3 + 12x2y + 5y5 - 10x2y3 - 15x4y
				break;	
		case 5: Z_mn = (   Poly(5,0,x,y)-10*Poly(3,2,x,y)+5*Poly(1,4,x,y));//Z(5,5) = y5 - 10x2y3 + 5x4y
				break;	
		}
		break;
	case 6:	 //Z(6,*)
	switch(n){
		case 0: Z_mn = (6*Poly(1,5,x,y)-20*Poly(3,3,x,y)+6*Poly(5,1,x,y));//Z(6,0) = 6x5y - 20x3y3 + 6xy5
				break;
		case 1: Z_mn = (20*Poly(1,3,x,y)-20*Poly(3,1,x,y)-24*Poly(1,5,x,y)+24*Poly(5,1,x,y));//Z(6,1) = 20x3y - 20xy3 - 24x5y + 24xy5
				break;	
		case 2: Z_mn = (12*Poly(1,1,x,y)-40*Poly(3,1,x,y)-40*Poly(1,3,x,y)+30*Poly(1,5,x,y)+60*Poly(3,3,x,y)+30*Poly(5,1,x,y));//Z(6,2) = 12xy - 40x3y - 40xy3 + 30x5y + 60x3y3 - 30xy5
				break;	
		case 3: Z_mn = (-Poly(0,0,x,y)+12*Poly(0,2,x,y)+12*Poly(2,0,x,y)-30*Poly(0,4,x,y)-60*Poly(2,2,x,y)-30*Poly(4,0,x,y)+20*Poly(0,6,x,y)+60*Poly(2,4,x,y)+60*Poly(4,2,x,y)+20*Poly(6,0,x,y));//Z(6,3) = -1 + 12x2 + 12y2 - 30x4 - 60x2y2 - 30y4 + 20x6 + 60x4y2 + 60x2y4 + 20y6
				break;	
		case 4: Z_mn = (-6*Poly(0,2,x,y)+6*Poly(2,0,x,y)+20*Poly(0,4,x,y)-20*Poly(4,0,x,y)-15*Poly(0,6,x,y)-15*Poly(2,4,x,y)+15*Poly(4,2,x,y)+15*Poly(6,0,x,y)); //Z(6,4) = -6x2 + 6y2 + 20x4 - 20y4 - 15x6 - 15x4y2 + 15x2y4 + 15y6
				break;	
		case 5: Z_mn = (-5*Poly(0,4,x,y)+30*Poly(2,2,x,y)-5*Poly(4,0,x,y)+6*Poly(0,6,x,y)-30*Poly(2,4,x,y)-30*Poly(4,2,x,y)+6*Poly(6,0,x,y));//Z(6,5) = -5x4 + 30x2y2 - 5y4 + 6x6 - 30x4y2 - 30x2y4 + 6y6
				break;	
		case 6: Z_mn = (-Poly(0,6,x,y)+15*Poly(2,4,x,y)-15*Poly(4,2,x,y)+Poly(6,0,x,y));//Z(6,6) = -x6 + 15x4y2 - 15x2y4 + y6
				break;
	}
		break;
	case 7:	 //Z(7,*)
	switch(n){
		case 0: Z_mn = (-Poly(0,7,x,y)+21*Poly(2,5,x,y)-35*Poly(4,3,x,y)+7*Poly(6,1,x,y));//Z(7,0) = -x7 + 21x5y2 - 35x3y4 + 7xy6
				break;	
		case 1: Z_mn = (-6*Poly(0,5,x,y) + 60*Poly(2,3,x,y) - 30*Poly(4,1,x,y)+7*Poly(0,7,x,y)- 63*Poly(2,5,x,y)- 35*Poly(4,3,x,y)+35*Poly(6,1,x,y));//Z(7,1) = -6x5 + 60x3y2 - 30xy4 + 7x7 - 63x5y2 - 35x3y4 + 35xy6
				break;	
		case 2: Z_mn = (-10*Poly(0,3,x,y) + 30*Poly(2,1,x,y) + 30*Poly(0,5,x,y)-60*Poly(2,3,x,y) -90*Poly(4,1,x,y) - 21*Poly(0,7,x,y) + 21*Poly(2,5,x,y) + 105*Poly(4,3,x,y)+63*Poly(6,1,x,y));//Z(7,2) = –10x3 + 30xy2 + 30x5 – 60x3y2 – 90xy4 – 21x7 + 21x5y2 + 105x3y4 + 63xy6
				break;	
		case 3: Z_mn = (-4*Poly(0,1,x,y) + 30*Poly(0,3,x,y) + 30*Poly(2,1,x,y)-60*Poly(0,5,x,y) - 120*Poly(2,3,x,y) -60*Poly(4,1,x,y)+35*Poly(0,7,x,y) + 105*Poly(2,5,x,y)+105*Poly(4,3,x,y)+35*Poly(6,1,x,y));//Z(7,3) = –4x + 30x3 + 30xy2 – 60x5 – 120x3y2 – 60xy4 + 35x7 + 105x5y2 + 105x3y4 + 35xy6
				break;	
		case 4: Z_mn = (-4*Poly(1,0,x,y)+30*Poly(3,0,x,y)+30*Poly(1,2,x,y)-60*Poly(5,0,x,y)-120*Poly(3,2,x,y)-60*Poly(1,4,x,y)+35*Poly(7,0,x,y)+105*Poly(5,2,x,y)+105*Poly(3,4,x,y)+35*Poly(1,6,x,y));	//Z(7,4) = –4y + 30y3 + 30x2y – 60y5 – 120x2y3 – 60x4y + 35y7 + 105x2y5 + 105x4y3 + 35x6y
				break;	
		case 5: Z_mn = ( 10*Poly(3,0,x,y) - 30*Poly(1,2,x,y) -30*Poly(5,0,x,y) + 60*Poly(3,2,x,y) + 90*Poly(1,4,x,y) + 21*Poly(7,0,x,y) -21*Poly(5,2,x,y) -105*Poly(3,4,x,y) -63*Poly(1,6,x,y) );//Z(7,5) = 10y3 – 30x2y – 30y5 + 60x2y3 + 90x4y + 21y7 – 21x2y5 – 105x4y3 + 63x6y
				break;	
		case 6: Z_mn = (-6*Poly(5,0,x,y) +60*Poly(3,2,x,y)-30*Poly(1,4,x,y)+7*Poly(7,0,x,y)-63*Poly(5,2,x,y)-35*Poly(3,4,x,y)+35*Poly(1,6,x,y));	//Z(7,6) = –6y5 + 60x2y3 – 30x4y + 7y7 – 63x2y5 – 35x4y3 + 35x6y
				break;	
		case 7: Z_mn = (Poly(7,0,x,y)-21*Poly(5,2,x,y)+35*Poly(3,4,x,y)-7*Poly(1,6,x,y));//Z(7,7) = y7 – 21x2y5 + 35x4y3 – 7x6y
				break;
		}
	break;
	case 8:	 //Z(8,*)
	switch(n){
		case 0: Z_mn = (-8*Poly(1,7,x,y)+56*Poly(3,5,x,y)-56*Poly(5,3,x,y)+8*Poly(7,1,x,y));//Z(8,0) = –8x7y + 56x5y3 – 56x3y5 + 8xy7
				break;	
		case 1: Z_mn = (-42*Poly(1,5,x,y)+140*Poly(3,3,x,y)-42*Poly(5,1,x,y)+48*Poly(1,7,x,y)-112*Poly(3,5,x,y)-112*Poly(5,3,x,y)+48*Poly(7,1,x,y));//Z(8,1) = –42x5y + 140x3y3 – 42xy5 + 48x7y – 112x5y3 – 112x3y5 + 48xy7
				break;	
		case 2: Z_mn = (-60*Poly(1,3,x,y)+60*Poly(3,1,x,y)+168*Poly(1,5,x,y)-168*Poly(5,1,x,y)-112*Poly(1,7,x,y)-112*Poly(3,5,x,y)+112*Poly(5,3,x,y)+112*Poly(7,1,x,y));//Z(8,2) = –60x3y + 60xy3 + 168x5y – 168xy5 – 112x7y – 112x5y3 + 112x3y5 + 112xy7
				break;	
		case 3: Z_mn = (-20*Poly(1,1,x,y)+120*Poly(1,3,x,y)+120*Poly(3,1,x,y)-210*Poly(1,5,x,y)-420*Poly(3,3,x,y)-210*Poly(5,1,x,y)+112*Poly(1,7,x,y)+336*Poly(3,5,x,y)+336*Poly(5,3,x,y)+112*Poly(7,1,x,y));//Z(8,3) = –20xy + 120x3y + 120xy3 – 210x5y – 420x3y3 – 210xy5 + 112x7y + 336x5y3 + 336x3y5 + 112xy7
				break;	
		case 4: Z_mn = (Poly(0,0,x,y)-20*Poly(0,2,x,y)-20*Poly(2,0,x,y)+90*Poly(0,4,x,y)+180*Poly(2,2,x,y)+90*Poly(4,0,x,y)-140*Poly(0,6,x,y)-420*Poly(2,4,x,y)-420*Poly(4,2,x,y)-140*Poly(6,0,x,y)+70*Poly(8,0,x,y)+280*Poly(2,6,x,y)+420*Poly(4,4,x,y)+280*Poly(6,2,x,y)+70*Poly(0,8,x,y));//Z(8,4) = 1 – 20x2 – 20y2 + 90x4 + 180x2y2 + 90y4 – 140x6 – 420x4y2 – 420x2y4 – 140y6 + 70x8 + 280x6y2 + 420x4y4 + 280x2y6 + 70y8
				break;	
		case 5: Z_mn = (10*Poly(0,2,x,y)-10*Poly(2,0,x,y)-60*Poly(0,4,x,y)+105*Poly(2,4,x,y)-105*Poly(4,2,x,y)+60*Poly(4,0,x,y)+105*Poly(0,6,x,y)-105*Poly(6,0,x,y)-56*Poly(0,8,x,y)-112*Poly(2,6,x,y)+112*Poly(6,2,x,y)+56*Poly(8,0,x,y));	//Z(8,5) = 10x2 – 10y2 – 60x4 + 105x4y2 – 105x2y4 + 60y4 + 105x6 – 105y6 – 56x8 – 112x6y2 + 112x2y6 + 56y8
				break;	
		case 6: Z_mn = (15*Poly(0,4,x,y)-90*Poly(2,2,x,y)+15*Poly(4,0,x,y)-42*Poly(0,6,x,y)+210*Poly(2,4,x,y)+210*Poly(4,2,x,y)-42*Poly(6,0,x,y)+28*Poly(0,8,x,y)-112*Poly(2,6,x,y)-280*Poly(4,4,x,y)-112*Poly(6,2,x,y)+28*Poly(8,0,x,y));//Z(8,6) = 15x4 – 90x2y2 + 15y4 – 42x6 + 210x4y2 + 210x2y4 – 42y6 + 28x8 – 112x6y2 – 280x4y4 – 112x2y6 + 28y8
				break;	
		case 7: Z_mn = (7*Poly(0,6,x,y)-105*Poly(2,4,x,y)+105*Poly(4,2,x,y)-7*Poly(6,0,x,y)-8*Poly(0,8,x,y)+112*Poly(2,6,x,y)-112*Poly(6,2,x,y)+8*Poly(8,0,x,y));//Z(8,7) = 7x6 – 105x4y2 + 105x2y4 – 7y6 – 8x8 + 112x6y2 – 112x2y6 + 8y8
				break;
		case 8: Z_mn = (Poly(0,8,x,y)-28*Poly(2,6,x,y)+70*Poly(4,4,x,y)-28*Poly(6,2,x,y)+1*Poly(8,0,x,y));//Z(8,8) = x8 – 28x6y2 + 70x4y4 – 28x2y6 + y8
				break;
		}
		break;
	default:
		Z_mn = 0.0;
		break;	
	}
	
	// add scalor value
	Z_mn = m_K[m][n]*Z_mn;

	return (Z_mn);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// This is the Partial Derivative of Zernike Polynomila of Z(m, n) with respect to x
//
////////////////////////////////////////////////////////////////////////////////////////////////
double LLWzernike::ZernikePoly_Dx(unsigned int m, unsigned int n, double x, double y)
{	
	if((m>ZERNIKE_ORDER) || (n>ZERNIKE_ORDER))
	{
		return 0.0;
	}

	double Z_Dx = 0.0;
	int i=0, j=0;

	switch(m)
	{
	case 0:	 //Z(0,*)
			Z_Dx = Poly_Dx(0,0,x, y);
		break;	
	case 1:	 //Z(1,*)
		switch(n){
		case 0:	Z_Dx = Poly_Dx(0,1, x, y);
			break;	
		case 1: Z_Dx = Poly_Dx(1,0, x, y);
			break;
		}
		break;	
	case 2:	 //Z(2,*)
		switch(n){
		case 0: Z_Dx = (2*Poly_Dx(1,1,x,y));	//Z(2,0) = 2xy
				break;	
		case 1: Z_Dx = (-Poly_Dx(0,0,x,y)+2*Poly_Dx(0,2,x,y)+2*Poly_Dx(2,0,x,y));	//Z(2,1) = 1+2x2 + 2y2
				break;	
		case 2: Z_Dx = (-Poly_Dx(0,2,x,y)+Poly_Dx(2,0,x,y));	//Z(2,2) = -x2 + y2
			break;
		}
		break;
	case 3:	 //Z(3,*)
		switch(n){
		case 0: Z_Dx = (-Poly_Dx(0,3,x,y)+3*Poly_Dx(2,1,x,y));// Z(3,0) = -x3 + 3xy2
				break;	
		case 1: Z_Dx = (-2*Poly_Dx(0,1,x,y)+3*Poly_Dx(0,3,x,y)+3*Poly_Dx(2,1,x,y));	// Z(3,1) = -2x + 3x3 + 3xy2 ((0),(-2,0),(0,0,0),(3,0,3,0))
				break;	
		case 2: Z_Dx = (-2*Poly_Dx(1,0,x,y)+3*Poly_Dx(3,0,x,y)+3*Poly_Dx(1,2,x,y));// Z(3,2) = -2y + 3y3 + 3x2y ((0),(0,-2),(0,0,0),(0,3,0,3))
				break;	
		case 3: Z_Dx = (Poly_Dx(3,0,x,y)-3*Poly_Dx(1,2,x,y));// Z(3,3) = y3 - 3x2y ((0),(0,0),(0,0,0),(0,-3,0,1))
				break;
		}
		break;
	case 4:	 //Z(4,*)
		switch(n){
		case 0: Z_Dx = (-4*Poly_Dx(1,3,x,y)+4*Poly_Dx(3,1,x,y));// Z(4,0) = -4x3y + 4xy3 ((0),(0,0),(0,0,0),(0,0,0,0),(0,-4,0,4,0))
				break;	
		case 1: Z_Dx = (-6*Poly_Dx(1,1,x,y)+8*Poly_Dx(1,3,x,y)+8*Poly_Dx(3,1,x,y));// Z(4,1) = -6xy + 8x3y + 8xy3 ((0),(0,0),(0,-6,0),(0,0,0,0),(0,8,0,8,0))
				break;	
		case 2: Z_Dx = (Poly_Dx(0,0,x,y)-6*Poly_Dx(0,2,x,y)-6*Poly_Dx(2,0,x,y)+6*Poly_Dx(0,4,x,y)+12*Poly_Dx(2,2,x,y)+6*Poly_Dx(4,0,x,y));// Z(4,2) = 1 - 6x2 - 6y2 + 6x4 + 12x2y2 + 6y4 ((1),(0,0),(-6,0,-6),(0,0,0,0),(6,0,12,0,6))
				break;	
		case 3: Z_Dx = (3*Poly_Dx(0,2,x,y)-3*Poly_Dx(2,0,x,y)-4*Poly_Dx(0,4,x,y)+4*Poly_Dx(4,0,x,y));// Z(4,3) = 3x2 - 3y2 - 4x4 + 4y4 ((0),(0,0),(3,0,-3),(0,0,0,0),(-4,0,0,0,4))
				break;	
		case 4: Z_Dx = (Poly_Dx(0,4,x,y)-6*Poly_Dx(2,2,x,y)+Poly_Dx(4,0,x,y));// Z(4,4) = x4 - 6x2y2 + y4 ((0),(0,0),(0,0,0),(0,0,0,0),(1,0,-6,0,1))
				break;
		}
		break;
	case 5:	 //Z(5,*)
		switch(n){
		case 0: Z_Dx = (   Poly_Dx(0,5,x,y)-10*Poly_Dx(2,3,x,y)+5*Poly_Dx(4,1,x,y));// Z(5,0) = x5 - 10x3y2 + 5xy4
				break;	
		case 1: Z_Dx = ( 4*Poly_Dx(0,3,x,y)-12*Poly_Dx(2,1,x,y)-5*Poly_Dx(0,5,x,y)+10*Poly_Dx(2,3,x,y)+15*Poly_Dx(4,1,x,y));//Z(5,1) = 4x3 - 12xy2 - 5x5 + 10x3y2 + 15xy4
				break;
		case 2: Z_Dx = ( 3*Poly_Dx(0,1,x,y)-12*Poly_Dx(0,3,x,y)-12*Poly_Dx(2,1,x,y)+10*Poly_Dx(0,5,x,y) + 20*Poly_Dx(2,3,x,y) + 10*Poly_Dx(4,1,x,y));//Z(5,2) = 3x - 12x3 - 12xy2 + 10x5 + 20x3y2 + 10xy4
				break;	
		case 3: Z_Dx = ( 3*Poly_Dx(1,0,x,y)-12*Poly_Dx(3,0,x,y)-12*Poly_Dx(1,2,x,y)+10*Poly_Dx(5,0,x,y)+20*Poly_Dx(3,2,x,y)+10*Poly_Dx(1,4,x,y));//Z(5,3) = 3y - 12y3 - 12x2y + 10y5 + 20x2y3 + 10x4y
				break;	
		case 4: Z_Dx = (-4*Poly_Dx(3,0,x,y)+12*Poly_Dx(1,2,x,y)+5*Poly_Dx(5,0,x,y)-10*Poly_Dx(3,2,x,y)-15*Poly_Dx(1,4,x,y));//Z(5,4) = -4y3 + 12x2y + 5y5 - 10x2y3 - 15x4y
				break;	
		case 5: Z_Dx = (   Poly_Dx(5,0,x,y)-10*Poly_Dx(3,2,x,y)+5*Poly_Dx(1,4,x,y));//Z(5,5) = y5 - 10x2y3 + 5x4y
				break;	
		}
		break;
	case 6:	 //Z(6,*)
	switch(n){
		case 0: Z_Dx = (6*Poly_Dx(1,5,x,y)-20*Poly_Dx(3,3,x,y)+6*Poly_Dx(5,1,x,y));//Z(6,0) = 6x5y - 20x3y3 + 6xy5
				break;
		case 1: Z_Dx = (20*Poly_Dx(1,3,x,y)-20*Poly_Dx(3,1,x,y)-24*Poly_Dx(1,5,x,y)+24*Poly_Dx(5,1,x,y));//Z(6,1) = 20x3y - 20xy3 - 24x5y + 24xy5
				break;	
		case 2: Z_Dx = (12*Poly_Dx(1,1,x,y)-40*Poly_Dx(3,1,x,y)-40*Poly_Dx(1,3,x,y)+30*Poly_Dx(1,5,x,y)+60*Poly_Dx(3,3,x,y)+30*Poly_Dx(5,1,x,y));//Z(6,2) = 12xy - 40x3y - 40xy3 + 30x5y + 60x3y3 - 30xy5
				break;	
		case 3: Z_Dx = (-Poly_Dx(0,0,x,y)+12*Poly_Dx(0,2,x,y)+12*Poly_Dx(2,0,x,y)-30*Poly_Dx(0,4,x,y)-60*Poly_Dx(2,2,x,y)-30*Poly_Dx(4,0,x,y)+20*Poly_Dx(0,6,x,y)+60*Poly_Dx(2,4,x,y)+60*Poly_Dx(4,2,x,y)+20*Poly_Dx(6,0,x,y));//Z(6,3) = -1 + 12x2 + 12y2 - 30x4 - 60x2y2 - 30y4 + 20x6 + 60x4y2 + 60x2y4 + 20y6
				break;	
		case 4: Z_Dx = (-6*Poly_Dx(0,2,x,y)+6*Poly_Dx(2,0,x,y)+20*Poly_Dx(0,4,x,y)-20*Poly_Dx(4,0,x,y)-15*Poly_Dx(0,6,x,y)-15*Poly_Dx(2,4,x,y)+15*Poly_Dx(4,2,x,y)+15*Poly_Dx(6,0,x,y)); //Z(6,4) = -6x2 + 6y2 + 20x4 - 20y4 - 15x6 - 15x4y2 + 15x2y4 + 15y6
				break;	
		case 5: Z_Dx = (-5*Poly_Dx(0,4,x,y)+30*Poly_Dx(2,2,x,y)-5*Poly_Dx(4,0,x,y)+6*Poly_Dx(0,6,x,y)-30*Poly_Dx(2,4,x,y)-30*Poly_Dx(4,2,x,y)+6*Poly_Dx(6,0,x,y));//Z(6,5) = -5x4 + 30x2y2 - 5y4 + 6x6 - 30x4y2 - 30x2y4 + 6y6
				break;	
		case 6: Z_Dx = (-Poly_Dx(0,6,x,y)+15*Poly_Dx(2,4,x,y)-15*Poly_Dx(4,2,x,y)+Poly_Dx(6,0,x,y));//Z(6,6) = -x6 + 15x4y2 - 15x2y4 + y6
				break;
	}
		break;
	case 7:	 //Z(7,*)
	switch(n){
		case 0: Z_Dx = (-Poly_Dx(0,7,x,y)+21*Poly_Dx(2,5,x,y)-35*Poly_Dx(4,3,x,y)+7*Poly_Dx(6,1,x,y));//Z(7,0) = -x7 + 21x5y2 - 35x3y4 + 7xy6
				break;	
		case 1: Z_Dx = (-6*Poly_Dx(0,5,x,y) + 60*Poly_Dx(2,3,x,y) - 30*Poly_Dx(4,1,x,y)+7*Poly_Dx(0,7,x,y)- 63*Poly_Dx(2,5,x,y)- 35*Poly_Dx(4,3,x,y)+35*Poly_Dx(6,1,x,y));//Z(7,1) = -6x5 + 60x3y2 - 30xy4 + 7x7 - 63x5y2 - 35x3y4 + 35xy6
				break;	
		case 2: Z_Dx = (-10*Poly_Dx(0,3,x,y) + 30*Poly_Dx(2,1,x,y) + 30*Poly_Dx(0,5,x,y)-60*Poly_Dx(2,3,x,y) -90*Poly_Dx(4,1,x,y) - 21*Poly_Dx(0,7,x,y) + 21*Poly_Dx(2,5,x,y) + 105*Poly_Dx(4,3,x,y)+63*Poly_Dx(6,1,x,y));//Z(7,2) = –10x3 + 30xy2 + 30x5 – 60x3y2 – 90xy4 – 21x7 + 21x5y2 + 105x3y4 + 63xy6
				break;	
		case 3: Z_Dx = (-4*Poly_Dx(0,1,x,y) + 30*Poly_Dx(0,3,x,y) + 30*Poly_Dx(2,1,x,y)-60*Poly_Dx(0,5,x,y) - 120*Poly_Dx(2,3,x,y) -60*Poly_Dx(4,1,x,y)+35*Poly_Dx(0,7,x,y) + 105*Poly_Dx(2,5,x,y)+105*Poly_Dx(4,3,x,y)+35*Poly_Dx(6,1,x,y));//Z(7,3) = –4x + 30x3 + 30xy2 – 60x5 – 120x3y2 – 60xy4 + 35x7 + 105x5y2 + 105x3y4 + 35xy6
				break;	
		case 4: Z_Dx = (-4*Poly_Dx(1,0,x,y)+30*Poly_Dx(3,0,x,y)+30*Poly_Dx(1,2,x,y)-60*Poly_Dx(5,0,x,y)-120*Poly_Dx(3,2,x,y)-60*Poly_Dx(1,4,x,y)+35*Poly_Dx(7,0,x,y)+105*Poly_Dx(5,2,x,y)+105*Poly_Dx(3,4,x,y)+35*Poly_Dx(1,6,x,y));	//Z(7,4) = –4y + 30y3 + 30x2y – 60y5 – 120x2y3 – 60x4y + 35y7 + 105x2y5 + 105x4y3 + 35x6y
				break;	
		case 5: Z_Dx = ( 10*Poly_Dx(3,0,x,y) - 30*Poly_Dx(1,2,x,y) -30*Poly_Dx(5,0,x,y) + 60*Poly_Dx(3,2,x,y) + 90*Poly_Dx(1,4,x,y) + 21*Poly_Dx(7,0,x,y) -21*Poly_Dx(5,2,x,y) -105*Poly_Dx(3,4,x,y) -63*Poly_Dx(1,6,x,y) );//Z(7,5) = 10y3 – 30x2y – 30y5 + 60x2y3 + 90x4y + 21y7 – 21x2y5 – 105x4y3 + 63x6y
				break;	
		case 6: Z_Dx = (-6*Poly_Dx(5,0,x,y) +60*Poly_Dx(3,2,x,y)-30*Poly_Dx(1,4,x,y)+7*Poly_Dx(7,0,x,y)-63*Poly_Dx(5,2,x,y)-35*Poly_Dx(3,4,x,y)+35*Poly_Dx(1,6,x,y));	//Z(7,6) = –6y5 + 60x2y3 – 30x4y + 7y7 – 63x2y5 – 35x4y3 + 35x6y
				break;	
		case 7: Z_Dx = (Poly_Dx(7,0,x,y)-21*Poly_Dx(5,2,x,y)+35*Poly_Dx(3,4,x,y)-7*Poly_Dx(1,6,x,y));//Z(7,7) = y7 – 21x2y5 + 35x4y3 – 7x6y
				break;
		}
	break;
	case 8:	 //Z(8,*)
	switch(n){
		case 0: Z_Dx = (-8*Poly_Dx(1,7,x,y)+56*Poly_Dx(3,5,x,y)-56*Poly_Dx(5,3,x,y)+8*Poly_Dx(7,1,x,y));//Z(8,0) = –8x7y + 56x5y3 – 56x3y5 + 8xy7
				break;	
		case 1: Z_Dx = (-42*Poly_Dx(1,5,x,y)+140*Poly_Dx(3,3,x,y)-42*Poly_Dx(5,1,x,y)+48*Poly_Dx(1,7,x,y)-112*Poly_Dx(3,5,x,y)-112*Poly_Dx(5,3,x,y)+48*Poly_Dx(7,1,x,y));//Z(8,1) = –42x5y + 140x3y3 – 42xy5 + 48x7y – 112x5y3 – 112x3y5 + 48xy7
				break;	
		case 2: Z_Dx = (-60*Poly_Dx(1,3,x,y)+60*Poly_Dx(3,1,x,y)+168*Poly_Dx(1,5,x,y)-168*Poly_Dx(5,1,x,y)-112*Poly_Dx(1,7,x,y)-112*Poly_Dx(3,5,x,y)+112*Poly_Dx(5,3,x,y)+112*Poly_Dx(7,1,x,y));//Z(8,2) = –60x3y + 60xy3 + 168x5y – 168xy5 – 112x7y – 112x5y3 + 112x3y5 + 112xy7
				break;	
		case 3: Z_Dx = (-20*Poly_Dx(1,1,x,y)+120*Poly_Dx(1,3,x,y)+120*Poly_Dx(3,1,x,y)-210*Poly_Dx(1,5,x,y)-420*Poly_Dx(3,3,x,y)-210*Poly_Dx(5,1,x,y)+112*Poly_Dx(1,7,x,y)+336*Poly_Dx(3,5,x,y)+336*Poly_Dx(5,3,x,y)+112*Poly_Dx(7,1,x,y));//Z(8,3) = –20xy + 120x3y + 120xy3 – 210x5y – 420x3y3 – 210xy5 + 112x7y + 336x5y3 + 336x3y5 + 112xy7
				break;	
		case 4: Z_Dx = (Poly_Dx(0,0,x,y)-20*Poly_Dx(0,2,x,y)-20*Poly_Dx(2,0,x,y)+90*Poly_Dx(0,4,x,y)+180*Poly_Dx(2,2,x,y)+90*Poly_Dx(4,0,x,y)-140*Poly_Dx(0,6,x,y)-420*Poly_Dx(2,4,x,y)-420*Poly_Dx(4,2,x,y)-140*Poly_Dx(6,0,x,y)+70*Poly_Dx(8,0,x,y)+280*Poly_Dx(2,6,x,y)+420*Poly_Dx(4,4,x,y)+280*Poly_Dx(6,2,x,y)+70*Poly_Dx(0,8,x,y));//Z(8,4) = 1 – 20x2 – 20y2 + 90x4 + 180x2y2 + 90y4 – 140x6 – 420x4y2 – 420x2y4 – 140y6 + 70x8 + 280x6y2 + 420x4y4 + 280x2y6 + 70y8
				break;	
		case 5: Z_Dx = (10*Poly_Dx(0,2,x,y)-10*Poly_Dx(2,0,x,y)-60*Poly_Dx(0,4,x,y)+105*Poly_Dx(2,4,x,y)-105*Poly_Dx(4,2,x,y)+60*Poly_Dx(4,0,x,y)+105*Poly_Dx(0,6,x,y)-105*Poly_Dx(6,0,x,y)-56*Poly_Dx(0,8,x,y)-112*Poly_Dx(2,6,x,y)+112*Poly_Dx(6,2,x,y)+56*Poly_Dx(8,0,x,y));	//Z(8,5) = 10x2 – 10y2 – 60x4 + 105x4y2 – 105x2y4 + 60y4 + 105x6 – 105y6 – 56x8 – 112x6y2 + 112x2y6 + 56y8
				break;	
		case 6: Z_Dx = (15*Poly_Dx(0,4,x,y)-90*Poly_Dx(2,2,x,y)+15*Poly_Dx(4,0,x,y)-42*Poly_Dx(0,6,x,y)+210*Poly_Dx(2,4,x,y)+210*Poly_Dx(4,2,x,y)-42*Poly_Dx(6,0,x,y)+28*Poly_Dx(0,8,x,y)-112*Poly_Dx(2,6,x,y)-280*Poly_Dx(4,4,x,y)-112*Poly_Dx(6,2,x,y)+28*Poly_Dx(8,0,x,y));//Z(8,6) = 15x4 – 90x2y2 + 15y4 – 42x6 + 210x4y2 + 210x2y4 – 42y6 + 28x8 – 112x6y2 – 280x4y4 – 112x2y6 + 28y8
				break;	
		case 7: Z_Dx = (7*Poly_Dx(0,6,x,y)-105*Poly_Dx(2,4,x,y)+105*Poly_Dx(4,2,x,y)-7*Poly_Dx(6,0,x,y)-8*Poly_Dx(0,8,x,y)+112*Poly_Dx(2,6,x,y)-112*Poly_Dx(6,2,x,y)+8*Poly_Dx(8,0,x,y));//Z(8,7) = 7x6 – 105x4y2 + 105x2y4 – 7y6 – 8x8 + 112x6y2 – 112x2y6 + 8y8
				break;
		case 8: Z_Dx = (Poly_Dx(0,8,x,y)-28*Poly_Dx(2,6,x,y)+70*Poly_Dx(4,4,x,y)-28*Poly_Dx(6,2,x,y)+1*Poly_Dx(8,0,x,y));//Z(8,8) = x8 – 28x6y2 + 70x4y4 – 28x2y6 + y8
				break;
		}
		break;
	default:
		Z_Dx = 0.0;
		break;	
	}
	
	// add scalor value
	Z_Dx = m_K[m][n]*Z_Dx;

	return (Z_Dx);
}


////////////////////////////////////////////////////////////////////////////////////////
//  This is the Partial Derivative of Zernike Polynomila of Z(m, n) with respect to y
//
//////////////////////////////////////////////////////////////////////////////////////////
double LLWzernike::ZernikePoly_Dy(unsigned int m, unsigned int n, double x, double y)
{
	if((m>ZERNIKE_ORDER) || (n>ZERNIKE_ORDER))
	{
		return 0.0;
	}

	double Z_Dy = 0.0;
	int i=0, j=0;

	switch(m)
	{
	case 0:	 //Z(0,*)
			Z_Dy = Poly_Dy(0,0,x, y);
		break;	
	case 1:	 //Z(1,*)
		switch(n){
		case 0:	Z_Dy = Poly_Dy(0,1, x, y);
			break;	
		case 1: Z_Dy = Poly_Dy(1,0, x, y);
			break;
		}
		break;	
	case 2:	 //Z(2,*)
		switch(n){
		case 0: Z_Dy = (2*Poly_Dy(1,1,x,y));	//Z(2,0) = 2xy
				break;	
		case 1: Z_Dy = (-Poly_Dy(0,0,x,y)+2*Poly_Dy(0,2,x,y)+2*Poly_Dy(2,0,x,y));	//Z(2,1) = 1+2x2 + 2y2
				break;	
		case 2: Z_Dy = (-Poly_Dy(0,2,x,y)+Poly_Dy(2,0,x,y));	//Z(2,2) = -x2 + y2
			break;
		}
		break;
	case 3:	 //Z(3,*)
		switch(n){
		case 0: Z_Dy = (-Poly_Dy(0,3,x,y)+3*Poly_Dy(2,1,x,y));// Z(3,0) = -x3 + 3xy2
				break;	
		case 1: Z_Dy = (-2*Poly_Dy(0,1,x,y)+3*Poly_Dy(0,3,x,y)+3*Poly_Dy(2,1,x,y));	// Z(3,1) = -2x + 3x3 + 3xy2 ((0),(-2,0),(0,0,0),(3,0,3,0))
				break;	
		case 2: Z_Dy = (-2*Poly_Dy(1,0,x,y)+3*Poly_Dy(3,0,x,y)+3*Poly_Dy(1,2,x,y));// Z(3,2) = -2y + 3y3 + 3x2y ((0),(0,-2),(0,0,0),(0,3,0,3))
				break;	
		case 3: Z_Dy = (Poly_Dy(3,0,x,y)-3*Poly_Dy(1,2,x,y));// Z(3,3) = y3 - 3x2y ((0),(0,0),(0,0,0),(0,-3,0,1))
				break;
		}
		break;
	case 4:	 //Z(4,*)
		switch(n){
		case 0: Z_Dy = (-4*Poly_Dy(1,3,x,y)+4*Poly_Dy(3,1,x,y));// Z(4,0) = -4x3y + 4xy3 ((0),(0,0),(0,0,0),(0,0,0,0),(0,-4,0,4,0))
				break;	
		case 1: Z_Dy = (-6*Poly_Dy(1,1,x,y)+8*Poly_Dy(1,3,x,y)+8*Poly_Dy(3,1,x,y));// Z(4,1) = -6xy + 8x3y + 8xy3 ((0),(0,0),(0,-6,0),(0,0,0,0),(0,8,0,8,0))
				break;	
		case 2: Z_Dy = (Poly_Dy(0,0,x,y)-6*Poly_Dy(0,2,x,y)-6*Poly_Dy(2,0,x,y)+6*Poly_Dy(0,4,x,y)+12*Poly_Dy(2,2,x,y)+6*Poly_Dy(4,0,x,y));// Z(4,2) = 1 - 6x2 - 6y2 + 6x4 + 12x2y2 + 6y4 ((1),(0,0),(-6,0,-6),(0,0,0,0),(6,0,12,0,6))
				break;	
		case 3: Z_Dy = (3*Poly_Dy(0,2,x,y)-3*Poly_Dy(2,0,x,y)-4*Poly_Dy(0,4,x,y)+4*Poly_Dy(4,0,x,y));// Z(4,3) = 3x2 - 3y2 - 4x4 + 4y4 ((0),(0,0),(3,0,-3),(0,0,0,0),(-4,0,0,0,4))
				break;	
		case 4: Z_Dy = (Poly_Dy(0,4,x,y)-6*Poly_Dy(2,2,x,y)+Poly_Dy(4,0,x,y));// Z(4,4) = x4 - 6x2y2 + y4 ((0),(0,0),(0,0,0),(0,0,0,0),(1,0,-6,0,1))
				break;
		}
		break;
	case 5:	 //Z(5,*)
		switch(n){
		case 0: Z_Dy = (   Poly_Dy(0,5,x,y)-10*Poly_Dy(2,3,x,y)+5*Poly_Dy(4,1,x,y));// Z(5,0) = x5 - 10x3y2 + 5xy4
				break;	
		case 1: Z_Dy = ( 4*Poly_Dy(0,3,x,y)-12*Poly_Dy(2,1,x,y)-5*Poly_Dy(0,5,x,y)+10*Poly_Dy(2,3,x,y)+15*Poly_Dy(4,1,x,y));//Z(5,1) = 4x3 - 12xy2 - 5x5 + 10x3y2 + 15xy4
				break;
		case 2: Z_Dy = ( 3*Poly_Dy(0,1,x,y)-12*Poly_Dy(0,3,x,y)-12*Poly_Dy(2,1,x,y)+10*Poly_Dy(0,5,x,y) + 20*Poly_Dy(2,3,x,y) + 10*Poly_Dy(4,1,x,y));//Z(5,2) = 3x - 12x3 - 12xy2 + 10x5 + 20x3y2 + 10xy4
				break;	
		case 3: Z_Dy = ( 3*Poly_Dy(1,0,x,y)-12*Poly_Dy(3,0,x,y)-12*Poly_Dy(1,2,x,y)+10*Poly_Dy(5,0,x,y)+20*Poly_Dy(3,2,x,y)+10*Poly_Dy(1,4,x,y));//Z(5,3) = 3y - 12y3 - 12x2y + 10y5 + 20x2y3 + 10x4y
				break;	
		case 4: Z_Dy = (-4*Poly_Dy(3,0,x,y)+12*Poly_Dy(1,2,x,y)+5*Poly_Dy(5,0,x,y)-10*Poly_Dy(3,2,x,y)-15*Poly_Dy(1,4,x,y));//Z(5,4) = -4y3 + 12x2y + 5y5 - 10x2y3 - 15x4y
				break;	
		case 5: Z_Dy = (   Poly_Dy(5,0,x,y)-10*Poly_Dy(3,2,x,y)+5*Poly_Dy(1,4,x,y));//Z(5,5) = y5 - 10x2y3 + 5x4y
				break;	
		}
		break;
	case 6:	 //Z(6,*)
	switch(n){
		case 0: Z_Dy = (6*Poly_Dy(1,5,x,y)-20*Poly_Dy(3,3,x,y)+6*Poly_Dy(5,1,x,y));//Z(6,0) = 6x5y - 20x3y3 + 6xy5
				break;
		case 1: Z_Dy = (20*Poly_Dy(1,3,x,y)-20*Poly_Dy(3,1,x,y)-24*Poly_Dy(1,5,x,y)+24*Poly_Dy(5,1,x,y));//Z(6,1) = 20x3y - 20xy3 - 24x5y + 24xy5
				break;	
		case 2: Z_Dy = (12*Poly_Dy(1,1,x,y)-40*Poly_Dy(3,1,x,y)-40*Poly_Dy(1,3,x,y)+30*Poly_Dy(1,5,x,y)+60*Poly_Dy(3,3,x,y)+30*Poly_Dy(5,1,x,y));//Z(6,2) = 12xy - 40x3y - 40xy3 + 30x5y + 60x3y3 - 30xy5
				break;	
		case 3: Z_Dy = (-Poly_Dy(0,0,x,y)+12*Poly_Dy(0,2,x,y)+12*Poly_Dy(2,0,x,y)-30*Poly_Dy(0,4,x,y)-60*Poly_Dy(2,2,x,y)-30*Poly_Dy(4,0,x,y)+20*Poly_Dy(0,6,x,y)+60*Poly_Dy(2,4,x,y)+60*Poly_Dy(4,2,x,y)+20*Poly_Dy(6,0,x,y));//Z(6,3) = -1 + 12x2 + 12y2 - 30x4 - 60x2y2 - 30y4 + 20x6 + 60x4y2 + 60x2y4 + 20y6
				break;	
		case 4: Z_Dy = (-6*Poly_Dy(0,2,x,y)+6*Poly_Dy(2,0,x,y)+20*Poly_Dy(0,4,x,y)-20*Poly_Dy(4,0,x,y)-15*Poly_Dy(0,6,x,y)-15*Poly_Dy(2,4,x,y)+15*Poly_Dy(4,2,x,y)+15*Poly_Dy(6,0,x,y)); //Z(6,4) = -6x2 + 6y2 + 20x4 - 20y4 - 15x6 - 15x4y2 + 15x2y4 + 15y6
				break;	
		case 5: Z_Dy = (-5*Poly_Dy(0,4,x,y)+30*Poly_Dy(2,2,x,y)-5*Poly_Dy(4,0,x,y)+6*Poly_Dy(0,6,x,y)-30*Poly_Dy(2,4,x,y)-30*Poly_Dy(4,2,x,y)+6*Poly_Dy(6,0,x,y));//Z(6,5) = -5x4 + 30x2y2 - 5y4 + 6x6 - 30x4y2 - 30x2y4 + 6y6
				break;	
		case 6: Z_Dy = (-Poly_Dy(0,6,x,y)+15*Poly_Dy(2,4,x,y)-15*Poly_Dy(4,2,x,y)+Poly_Dy(6,0,x,y));//Z(6,6) = -x6 + 15x4y2 - 15x2y4 + y6
				break;
	}
		break;
	case 7:	 //Z(7,*)
	switch(n){
		case 0: Z_Dy = (-Poly_Dy(0,7,x,y)+21*Poly_Dy(2,5,x,y)-35*Poly_Dy(4,3,x,y)+7*Poly_Dy(6,1,x,y));//Z(7,0) = -x7 + 21x5y2 - 35x3y4 + 7xy6
				break;	
		case 1: Z_Dy = (-6*Poly_Dy(0,5,x,y) + 60*Poly_Dy(2,3,x,y) - 30*Poly_Dy(4,1,x,y)+7*Poly_Dy(0,7,x,y)- 63*Poly_Dy(2,5,x,y)- 35*Poly_Dy(4,3,x,y)+35*Poly_Dy(6,1,x,y));//Z(7,1) = -6x5 + 60x3y2 - 30xy4 + 7x7 - 63x5y2 - 35x3y4 + 35xy6
				break;	
		case 2: Z_Dy = (-10*Poly_Dy(0,3,x,y) + 30*Poly_Dy(2,1,x,y) + 30*Poly_Dy(0,5,x,y)-60*Poly_Dy(2,3,x,y) -90*Poly_Dy(4,1,x,y) - 21*Poly_Dy(0,7,x,y) + 21*Poly_Dy(2,5,x,y) + 105*Poly_Dy(4,3,x,y)+63*Poly_Dy(6,1,x,y));//Z(7,2) = –10x3 + 30xy2 + 30x5 – 60x3y2 – 90xy4 – 21x7 + 21x5y2 + 105x3y4 + 63xy6
				break;	
		case 3: Z_Dy = (-4*Poly_Dy(0,1,x,y) + 30*Poly_Dy(0,3,x,y) + 30*Poly_Dy(2,1,x,y)-60*Poly_Dy(0,5,x,y) - 120*Poly_Dy(2,3,x,y) -60*Poly_Dy(4,1,x,y)+35*Poly_Dy(0,7,x,y) + 105*Poly_Dy(2,5,x,y)+105*Poly_Dy(4,3,x,y)+35*Poly_Dy(6,1,x,y));//Z(7,3) = –4x + 30x3 + 30xy2 – 60x5 – 120x3y2 – 60xy4 + 35x7 + 105x5y2 + 105x3y4 + 35xy6
				break;	
		case 4: Z_Dy = (-4*Poly_Dy(1,0,x,y)+30*Poly_Dy(3,0,x,y)+30*Poly_Dy(1,2,x,y)-60*Poly_Dy(5,0,x,y)-120*Poly_Dy(3,2,x,y)-60*Poly_Dy(1,4,x,y)+35*Poly_Dy(7,0,x,y)+105*Poly_Dy(5,2,x,y)+105*Poly_Dy(3,4,x,y)+35*Poly_Dy(1,6,x,y));	//Z(7,4) = –4y + 30y3 + 30x2y – 60y5 – 120x2y3 – 60x4y + 35y7 + 105x2y5 + 105x4y3 + 35x6y
				break;	
		case 5: Z_Dy = ( 10*Poly_Dy(3,0,x,y) - 30*Poly_Dy(1,2,x,y) -30*Poly_Dy(5,0,x,y) + 60*Poly_Dy(3,2,x,y) + 90*Poly_Dy(1,4,x,y) + 21*Poly_Dy(7,0,x,y) -21*Poly_Dy(5,2,x,y) -105*Poly_Dy(3,4,x,y) -63*Poly_Dy(1,6,x,y) );//Z(7,5) = 10y3 – 30x2y – 30y5 + 60x2y3 + 90x4y + 21y7 – 21x2y5 – 105x4y3 + 63x6y
				break;	
		case 6: Z_Dy = (-6*Poly_Dy(5,0,x,y) +60*Poly_Dy(3,2,x,y)-30*Poly_Dy(1,4,x,y)+7*Poly_Dy(7,0,x,y)-63*Poly_Dy(5,2,x,y)-35*Poly_Dy(3,4,x,y)+35*Poly_Dy(1,6,x,y));	//Z(7,6) = –6y5 + 60x2y3 – 30x4y + 7y7 – 63x2y5 – 35x4y3 + 35x6y
				break;	
		case 7: Z_Dy = (Poly_Dy(7,0,x,y)-21*Poly_Dy(5,2,x,y)+35*Poly_Dy(3,4,x,y)-7*Poly_Dy(1,6,x,y));//Z(7,7) = y7 – 21x2y5 + 35x4y3 – 7x6y
				break;
		}
	break;
	case 8:	 //Z(8,*)
	switch(n){
		case 0: Z_Dy = (-8*Poly_Dy(1,7,x,y)+56*Poly_Dy(3,5,x,y)-56*Poly_Dy(5,3,x,y)+8*Poly_Dy(7,1,x,y));//Z(8,0) = –8x7y + 56x5y3 – 56x3y5 + 8xy7
				break;	
		case 1: Z_Dy = (-42*Poly_Dy(1,5,x,y)+140*Poly_Dy(3,3,x,y)-42*Poly_Dy(5,1,x,y)+48*Poly_Dy(1,7,x,y)-112*Poly_Dy(3,5,x,y)-112*Poly_Dy(5,3,x,y)+48*Poly_Dy(7,1,x,y));//Z(8,1) = –42x5y + 140x3y3 – 42xy5 + 48x7y – 112x5y3 – 112x3y5 + 48xy7
				break;	
		case 2: Z_Dy = (-60*Poly_Dy(1,3,x,y)+60*Poly_Dy(3,1,x,y)+168*Poly_Dy(1,5,x,y)-168*Poly_Dy(5,1,x,y)-112*Poly_Dy(1,7,x,y)-112*Poly_Dy(3,5,x,y)+112*Poly_Dy(5,3,x,y)+112*Poly_Dy(7,1,x,y));//Z(8,2) = –60x3y + 60xy3 + 168x5y – 168xy5 – 112x7y – 112x5y3 + 112x3y5 + 112xy7
				break;	
		case 3: Z_Dy = (-20*Poly_Dy(1,1,x,y)+120*Poly_Dy(1,3,x,y)+120*Poly_Dy(3,1,x,y)-210*Poly_Dy(1,5,x,y)-420*Poly_Dy(3,3,x,y)-210*Poly_Dy(5,1,x,y)+112*Poly_Dy(1,7,x,y)+336*Poly_Dy(3,5,x,y)+336*Poly_Dy(5,3,x,y)+112*Poly_Dy(7,1,x,y));//Z(8,3) = –20xy + 120x3y + 120xy3 – 210x5y – 420x3y3 – 210xy5 + 112x7y + 336x5y3 + 336x3y5 + 112xy7
				break;	
		case 4: Z_Dy = (Poly_Dy(0,0,x,y)-20*Poly_Dy(0,2,x,y)-20*Poly_Dy(2,0,x,y)+90*Poly_Dy(0,4,x,y)+180*Poly_Dy(2,2,x,y)+90*Poly_Dy(4,0,x,y)-140*Poly_Dy(0,6,x,y)-420*Poly_Dy(2,4,x,y)-420*Poly_Dy(4,2,x,y)-140*Poly_Dy(6,0,x,y)+70*Poly_Dy(8,0,x,y)+280*Poly_Dy(2,6,x,y)+420*Poly_Dy(4,4,x,y)+280*Poly_Dy(6,2,x,y)+70*Poly_Dy(0,8,x,y));//Z(8,4) = 1 – 20x2 – 20y2 + 90x4 + 180x2y2 + 90y4 – 140x6 – 420x4y2 – 420x2y4 – 140y6 + 70x8 + 280x6y2 + 420x4y4 + 280x2y6 + 70y8
				break;	
		case 5: Z_Dy = (10*Poly_Dy(0,2,x,y)-10*Poly_Dy(2,0,x,y)-60*Poly_Dy(0,4,x,y)+105*Poly_Dy(2,4,x,y)-105*Poly_Dy(4,2,x,y)+60*Poly_Dy(4,0,x,y)+105*Poly_Dy(0,6,x,y)-105*Poly_Dy(6,0,x,y)-56*Poly_Dy(0,8,x,y)-112*Poly_Dy(2,6,x,y)+112*Poly_Dy(6,2,x,y)+56*Poly_Dy(8,0,x,y));	//Z(8,5) = 10x2 – 10y2 – 60x4 + 105x4y2 – 105x2y4 + 60y4 + 105x6 – 105y6 – 56x8 – 112x6y2 + 112x2y6 + 56y8
				break;	
		case 6: Z_Dy = (15*Poly_Dy(0,4,x,y)-90*Poly_Dy(2,2,x,y)+15*Poly_Dy(4,0,x,y)-42*Poly_Dy(0,6,x,y)+210*Poly_Dy(2,4,x,y)+210*Poly_Dy(4,2,x,y)-42*Poly_Dy(6,0,x,y)+28*Poly_Dy(0,8,x,y)-112*Poly_Dy(2,6,x,y)-280*Poly_Dy(4,4,x,y)-112*Poly_Dy(6,2,x,y)+28*Poly_Dy(8,0,x,y));//Z(8,6) = 15x4 – 90x2y2 + 15y4 – 42x6 + 210x4y2 + 210x2y4 – 42y6 + 28x8 – 112x6y2 – 280x4y4 – 112x2y6 + 28y8
				break;	
		case 7: Z_Dy = (7*Poly_Dy(0,6,x,y)-105*Poly_Dy(2,4,x,y)+105*Poly_Dy(4,2,x,y)-7*Poly_Dy(6,0,x,y)-8*Poly_Dy(0,8,x,y)+112*Poly_Dy(2,6,x,y)-112*Poly_Dy(6,2,x,y)+8*Poly_Dy(8,0,x,y));//Z(8,7) = 7x6 – 105x4y2 + 105x2y4 – 7y6 – 8x8 + 112x6y2 – 112x2y6 + 8y8
				break;
		case 8: Z_Dy = (Poly_Dy(0,8,x,y)-28*Poly_Dy(2,6,x,y)+70*Poly_Dy(4,4,x,y)-28*Poly_Dy(6,2,x,y)+1*Poly_Dy(8,0,x,y));//Z(8,8) = x8 – 28x6y2 + 70x4y4 – 28x2y6 + y8
				break;
		}
		break;
	default:
		Z_Dy = 0.0;
		break;	
	}
	
	// add scalor value
	Z_Dy = m_K[m][n]*Z_Dy;

	return (Z_Dy);
}


//Zernike Polynonial expressions

/*
Z(0,0) = 1 
Z(1,0) = x
Z(1,1) = y
Z(2,0) = 2xy
Z(2,1) = 1+2x2 + 2y2
Z(2,2) = -x2 + y2
Z(3,0) = -x3 + 3xy2
Z(3,1) = -2x + 3x3 + 3xy2 ((0),(-2,0),(0,0,0),(3,0,3,0))
Z(3,2) = -2y + 3y3 + 3x2y ((0),(0,-2),(0,0,0),(0,3,0,3))
Z(3,3) = y3 - 3x2y ((0),(0,0),(0,0,0),(0,-3,0,1))
Z(4,0) = -4x3y + 4xy3 ((0),(0,0),(0,0,0),(0,0,0,0),(0,-4,0,4,0))
Z(4,1) = -6xy + 8x3y + 8xy3 ((0),(0,0),(0,-6,0),(0,0,0,0),(0,8,0,8,0))
Z(4,2) = 1 - 6x2 - 6y2 + 6x4 + 12x2y2 + 6y4 ((1),(0,0),(-6,0,-6),(0,0,0,0),(6,0,12,0,6))
Z(4,3) = 3x2 - 3y2 - 4x4 + 4y4 ((0),(0,0),(3,0,-3),(0,0,0,0),(-4,0,0,0,4))
Z(4,4) = x4 - 6x2y2 + y4 ((0),(0,0),(0,0,0),(0,0,0,0),(1,0,-6,0,1))
Z(5,0) = x5 - 10x3y2 + 5xy4
Z(5,1) = 4x3 - 12xy2 - 5x5 + 10x3y2 + 15xy4
Z(5,2) = 3x - 12x3 - 12xy2 + 10x5 + 20x3y2 + 10xy4
Z(5,3) = 3y - 12y3 - 12x2y + 10y5 + 20x2y3 + 10x4y
Z(5,4) = -4y3 + 12x2y + 5y5 - 10x2y3 - 15x4y
Z(5,5) = y5 - 10x2y3 + 5x4y
Z(6,0) = 6x5y - 20x3y3 + 6xy5
Z(6,1) = 20x3y - 20xy3 - 24x5y + 24xy5
Z(6,2) = 12xy - 40x3y - 40xy3 + 30x5y + 60x3y3 - 30xy5
Z(6,3) = -1 + 12x2 + 12y2 - 30x4 - 60x2y2 - 30y4 + 20x6 + 60x4y2 + 60x2y4 + 20y6
Z(6,4) = -6x2 + 6y2 + 20x4 - 20y4 - 15x6 - 15x4y2 + 15x2y4 + 15y6
Z(6,5) = -5x4 + 30x2y2 - 5y4 + 6x6 - 30x4y2 - 30x2y4 + 6y6
Z(6,6) = -x6 + 15x4y2 - 15x2y4 + y6
Z(7,0) = -x7 + 21x5y2 - 35x3y4 + 7xy6
Z(7,1) = -6x5 + 60x3y2 - 30xy4 + 7x7 - 63x5y2 - 35x3y4 + 35xy6
Z(7,2) = –10x3 + 30xy2 + 30x5 – 60x3y2 – 90xy4 – 21x7 + 21x5y2 + 105x3y4 + 63xy6
Z(7,3) = –4x + 30x3 + 30xy2 – 60x5 – 120x3y2 – 60xy4 + 35x7 + 105x5y2 + 105x3y4 + 35xy6
Z(7,4) = –4y + 30y3 + 30x2y – 60y5 – 120x2y3 – 60x4y + 35y7 + 105x2y5 + 105x4y3 + 35x6y
Z(7,5) = 10y3 – 30x2y – 30y5 + 60x2y3 + 90x4y + 21y7 – 21x2y5 – 105x4y3 + 63x6y
Z(7,6) = –6y5 + 60x2y3 – 30x4y + 7y7 – 63x2y5 – 35x4y3 + 35x6y
Z(7,7) = y7 – 21x2y5 + 35x4y3 – 7x6y
Z(8,0) = –8x7y + 56x5y3 – 56x3y5 + 8xy7
Z(8,1) = –42x5y + 140x3y3 – 42xy5 + 48x7y – 112x5y3 – 112x3y5 + 48xy7
Z(8,2) = –60x3y + 60xy3 + 168x5y – 168xy5 – 112x7y – 112x5y3 + 112x3y5 + 112xy7
Z(8,3) = –20xy + 120x3y + 120xy3 – 210x5y – 420x3y3 – 210xy5 – 112x7y + 336x5y3 + 336x3y5 + 112xy7
Z(8,4) = 1 – 20x2 – 20y2 + 90x4 + 180x2y2 + 90y4 – 140x6 – 420x4y2 – 420x2y4 – 140y6 + 70x8 + 280x6y2 + 420x4y4 + 280x2y6 + 70y8
Z(8,5) = 10x2 – 10y2 – 60x4 + 105x4y2 – 105x2y4 + 60y4 + 105x6 – 105y6 – 56x8 – 112x6y2 + 112x2y6 + 56y8
Z(8,6) = 15x4 – 90x2y2 + 15y4 – 42x6 + 210x4y2 + 210x2y4 – 42y6 + 28x8 – 112x6y2 – 280x4y4 – 112x2y6 + 28y8
Z(8,7) = 7x6 – 105x4y2 + 105x2y4 – 7y6 – 8x8 + 112x6y2 – 112x2y6 + 8y8
Z(8,8) = x8 – 28x6y2 + 70x4y4 – 28x2y6 + y8
Z(9,0) = x9 – 36x7y2 + 126x5y4 – 84x3y6 + 9xy8
Z(9,1) = 8x7 – 168x5y2 + 280x3y4 – 56xy6 – 9x9 + 180x7y2 – 126x5y4 – 252x3y6 + 63xy8
Z(9,2) = 21x5 – 210x3y2 + 105xy4 – 56x7 + 504x5y2 + 280x3y4 – 280xy6 + 36x9 – 288x7y2 – 504x5y4 + 180xy8
Z(9,3) = 20x3 – 60xy2 – 105x5 + 210x3y2 + 315xy4 + 168x7 – 168x5y2 – 840x3y4 – 504xy6 – 84x9 + 504x5y4 + 672x3y6 + 252xy8
Z(9,4) = 5x – 60x3 – 60xy2 + 210x5 + 420x3y2 + 210xy4 – 280x7 – 840x5y2 – 840x3y4 – 280xy6 + 126x9 + 504x7y2 + 756x5y4 + 504x3y6 + 126xy8
Z(9,5) = 5y – 60y3 – 60x2y + 210y5 + 420x2y3 + 210x4y – 280y7 – 840x2y5 – 840x4y3 – 280x6y + 126y9 + 504x2y7 + 756x4y5 + 504x6y3 + 126x8y
Z(9,6) = –20y3 + 60x2y + 105y5 – 210x2y3 – 315x4y – 168y7 + 168x2y5 + 840x4y3 + 504x6y + 84y9 – 504x4y5 – 672x6y3 – 252x8y
Z(9,7) = 21y5 – 210x2y3 + 105x4y – 56y7 + 504x2y5 + 280x4y3 – 280x6y + 36y9 – 288x2y7 – 504x4y5 + 180x8y
Z(9,8) = –8y7 + 168x2y5 – 280x4y3 + 56x6y + 9y9 – 180x2y7 + 126x4y5 + 252x6y3 – 63x8y
Z(9,9) = y9 – 36x2y7 + 126x4y5 – 84x6y3 + 9x8y
*/


