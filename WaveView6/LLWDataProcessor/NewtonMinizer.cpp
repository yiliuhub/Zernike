#include "stdafx.h"
#include <stdio.h>
#include <math.h>

#include "Marquard.h"

typedef double (*Function)(double*, double*, double*);


/*************************************************************************
*  int newton_minimize()	Newton's method to find minimum point of a fuction
*
*		fun_der1 --- the first order derivative function
*		fun_der2 --- the second order derivative function
*		*parameters --- the parameters used in the derivative functions to 
*						evaluate the derivatives
*		init_value --- the initial poisition to start search
*		*minimizer --- the minimum point
*		tolerence --- iteration tolerence
*
*  formula:
*												first_derivative
*			nest_position = current_position - -------------------
*												second_derivative
*
*  return value:
*			0 --- bad value, searching fails
*			1 --- successed, meets tolerence
*			2 --- meet maximum iteration number
*
*  Yi Liu  5/25/99
*
***************************************************************************/
int newton_minimize(double (*fun_der1)(double*, double),
					double (*fun_der2)(double*, double),
					double *parameters,
					double *minimizer,
					double init_value,
					double tolerence)
{
	int max_itr_num = 100;
	int rtn_value = 0;
	double delta;
	int i;
	double der1;
	double der2;
	double t;

	t = init_value;

	
	for(i = 0; i < max_itr_num; i++)
	{
		der1 = fun_der1(parameters, t);
		der2 = fun_der2(parameters, t);

		if(fabs(der2) <= EPS)
		{
			rtn_value = 0;
		}
		else
		{
			delta = der1/der2;
			t = t - delta;

			if(fabs(delta) < tolerence)
			{
				rtn_value = 1;
				break;
			}
		}
	}

	*minimizer = t;

	if(i >= max_itr_num)
		rtn_value = 2;


	return rtn_value;
}
