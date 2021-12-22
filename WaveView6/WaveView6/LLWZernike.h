/**************************************************************
*
*	Author:	Yi Liu
*	Date:	October 2010
*
****************************************************************/
#ifndef LLWZernike_H
#define LLWZernike_H

#include <stdio.h>
//#include <math.h>
//#include "ogpIpDefs.h"
//#include "OGPipCentroidDefs.h"
//#include "Geometry\geometry.h"

#define PIXEL_MICRON_RATIO		9
#define ZERNIKE_ORDER			8

// Hard coded coefficients for each single Zernike term
// !!!!!!!!!!!!   DON'T MAKE ANY CHANGE OF THE FOLOWING NUMBERS !!!!!!!!!!!!!!!!!!!!
/*
	m_K[0][0]=sqrt(2.);		m_K[0][1] = 0;			m_K[0][2] = 0;				m_K[0][3] = 0;			m_K[0][4] = 0;			m_K[0][5] = 0;			m_K[0][6] = 0;			m_K[0][7] = 0;			m_K[0][8] = 0;	
	m_K[1][0]=2;			m_K[1][1] = 2;			m_K[1][2] = 0;				m_K[1][3] = 0;			m_K[1][4] = 0;			m_K[1][5] = 0;			m_K[1][6] = 0;			m_K[1][7] = 0;			m_K[1][8] = 0;
	m_K[2][0]=sqrt(6.);		m_K[2][1]=sqrt(3.);		m_K[2][2]=sqrt(6.);			m_K[2][3] = 0;			m_K[2][4] = 0;			m_K[2][5] = 0;			m_K[2][6] = 0;			m_K[2][7] = 0;			m_K[2][8] = 0;
	m_K[3][0]=2*sqrt(2.);	m_K[3][1]=2*sqrt(2.);	m_K[3][2]=2*sqrt(2.);		m_K[3][3]=2*sqrt(2.);	m_K[3][4] = 0;			m_K[3][5] = 0;			m_K[3][6] = 0;			m_K[3][7] = 0;			m_K[3][8] = 0;
	m_K[4][0]=sqrt(10.);	m_K[4][1]=sqrt(10.);	m_K[4][2]=sqrt(5.);			m_K[4][3]=sqrt(10.);	m_K[4][4]=sqrt(10.);	m_K[4][5] = 0;			m_K[4][6] = 0;			m_K[4][7] = 0;			m_K[4][8] = 0;
	m_K[5][0]=2*sqrt(3.);	m_K[5][1]=2*sqrt(3.);	m_K[5][2]=2*sqrt(3.);		m_K[5][3]=2*sqrt(3.);	m_K[5][4]=2*sqrt(3.);	m_K[5][5]=2*sqrt(3.);	m_K[5][6] = 0;			m_K[5][7] = 0;			m_K[5][8] = 0;
	m_K[6][0]=sqrt(14.);	m_K[6][1]=sqrt(14.);	m_K[6][2]=sqrt(14.);		m_K[6][3]=sqrt(7.);		m_K[6][4]=sqrt(14.);	m_K[6][5]=sqrt(14.);	m_K[6][6]=sqrt(14.);	m_K[6][7] = 0;			m_K[6][8] = 0;
	m_K[7][0]=4;			m_K[7][1]=4;			m_K[7][2]=4;				m_K[7][3]=4;			m_K[7][4]=4;			m_K[7][5]=4;			m_K[7][6]=4;			m_K[7][7]=4;			m_K[7][8] = 0;
	m_K[8][0]=3*sqrt(2.);	m_K[8][1]=3*sqrt(2.);	m_K[8][2]=3*sqrt(2.);		m_K[8][3]=3*sqrt(2.);	m_K[8][4]=3;			m_K[8][5]=3*sqrt(2.);	m_K[8][6]=3*sqrt(2.);	m_K[8][7]=3*sqrt(2.);	m_K[8][8]=3*sqrt(2.);
*/
//////////////////////////////////////////////////////////////////////////
// this is hard coded coefficiient matrix of Zernike Polynomial
// Zernike(x, y) = sum(i,j)(m_K[i][j]*Zernike(i, j, x, y)
//////////////////////////////////////////////////////////////////////////
const double  m_K[ZERNIKE_ORDER+1][ZERNIKE_ORDER+1] =
{	
	{1.4142135623730951,	0.0000000000000000,	0.0000000000000000,		0.0000000000000000,		0.0000000000000000,		0.0000000000000000,		0.0000000000000000,	0.0000000000000000,	0.0000000000000000},
	{2.0000000000000000,	2.0000000000000000,	0.0000000000000000,		0.0000000000000000,		0.0000000000000000,		0.0000000000000000,		0.0000000000000000,	0.0000000000000000,	0.0000000000000000},
	{2.4494897427831779,	1.7320508075688772,	2.4494897427831779,		0.0000000000000000,		0.0000000000000000,		0.0000000000000000,		0.0000000000000000,	0.0000000000000000,	0.0000000000000000},
	{2*1.41421356237310,	2*1.41421356237310,	2*1.41421356237310,		2*1.41421356237310,		0.0000000000000000,		0.0000000000000000,		0.0000000000000000,	0.0000000000000000,	0.0000000000000000},
	{3.1622776601683795,	3.1622776601683795,	2.23606797749979,		3.1622776601683795,		3.1622776601683795,		0.0000000000000000,		0.0000000000000000,	0.0000000000000000,	0.0000000000000000},
	{2*1.73205080756888,	2*1.73205080756888,	2*1.7320508075688772,	2*1.73205080756888,		2*1.73205080756887,		2*1.73205080756888,		0.0000000000000000,	0.0000000000000000,	0.0000000000000000},
	{3.7416573867739413,	3.7416573867739413,	3.7416573867739413,		2.6457513110645907,		3.7416573867739413,		3.7416573867739413,		3.7416573867739413,	0.0000000000000000,	0.0000000000000000},
	{4.0000000000000000,	4.0000000000000000,	4.0000000000000000,		4.0000000000000000,		4.0000000000000000,		4.0000000000000000,		4.0000000000000000,	4.0000000000000000,	0.0000000000000000},
	{3*1.41421356237310,	3*1.41421356237310,	3*1.41421356237310,		3*1.41421356237310,		3.0000000000000000,		3*1.4142135623730951,	3*1.41421356237310,	3*1.41421356237310,	3*1.41421356237310}
};

//
//enum LLW_Error {
//	NO_ERROR = 0,
//	NO_ENOUGH_DATA,
//	INVALID_DATA,
//	INVALID_ORDER,
//	INVALID_FOCUS_LENGTH,
//	FAILED_OPEN_FILE,
//	MEMORY_ALLOCATION_FAILURE,
//	GENERAL_ERRO
//};


class LLWzernike
{
public:
	LLWzernike();
	~LLWzernike();
	int SetZernikeOrder(const unsigned int order);
	int GetZernikeOrder();
	int GetCoefficient(unsigned int m, unsigned int n, double& c); // get coefficient of Z(m,n)
	double Zernike(double x, double y);
	bool SetPointList(unsigned short size, double* pX, double* pY, double* pDraftX, double* pDraftY);
	bool SetFocusLength(double focus);
	bool SetOrder(int order);
	bool LoadDataFromFile();
	bool GoZernike();
	void Reset();

private:
	void TestDataIOFunction();
	double Poly(double order_x, double order_y, double x, double y);
	double Poly_Dx(double order_x, double order_y, double x, double y);
	double Poly_Dy(double order_x, double order_y, double x, double y);
	double ZernikePoly(unsigned int m, unsigned int n, double x, double y);
	double ZernikePoly_Dx(unsigned int m, unsigned int n, double x, double y);
	double ZernikePoly_Dy(unsigned int m, unsigned int n, double x, double y);
	double aPlusb(double a, double b);
	bool ValidityCheck();
	bool SetLinearSystem();
	bool SetNormalScalor(double ratio);
	bool PreProcessData();
	bool MapIJtoK(int i, int j, int& k);
	bool MapKtoIJ(int k, int& i, int& j);

private:	
	double m_dCoefficient[ZERNIKE_ORDER+1][ZERNIKE_ORDER+1];
	double** m_ppMatrix;
	double* m_pB;
	double* m_pX;
	double* m_pRefX;
	double* m_pRefY;
	double* m_pDraftX;
	double* m_pDraftY;
	double* m_pDeltaX;
	double* m_pDeltaY;
	double m_dCtrX;	// center X
	double m_dCtrY; // center Y
	int m_nSize;	// size of samples
	int m_nOrder;
	double m_dF;	// foucs length
	double m_dScalor; // normalization factor in pixel

//	Error_Type m_error;
};

#endif