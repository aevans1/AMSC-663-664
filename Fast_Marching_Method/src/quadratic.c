#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "fmm.h"

//solves a quadratic ax^2 + bx + c for real coefficients a,b,c, returns answer
//with list of REAL roots, returns NAN if complex
void solve_quadratic(double a, double b, double c, double* roots)
{
	//calculate discriminant
	double d;
	d = b*b - 4*a*c;

	if (d < 0)
	{
		//printf("Complex roots! Returning nAn \n");
	}
	//NOTE: d < 0, complex roots, will produce output NAN for both roots
	roots[0] = (-1*b - sqrt(d))/(2*a);
	roots[1] = (-1*b + sqrt(d))/(2*a);

	
}
/*Main for testing*/
//void main()
//{
//	double a;
//	double b;
//	double c;
//	double complex roots[2];
//
//	a = 1300.1119;
//	b = -69.1;
//	c = 2225.3;
//	solve_quadratic(a,b,c,roots);
//	printf("root 1 = %f + %f*I \n",creal(roots[0]),cimag(roots[0]));
//	printf("root 2 = %f + %f*I \n",creal(roots[1]),cimag(roots[1]));
//}
