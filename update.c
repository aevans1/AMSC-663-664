#include <stdio.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>
#include "fmm.h"
/*Find neighbors of point p in domain A
  Compute neighbors and eikonal update */

//input point p to update,
//array of poitns A(Ny is a global macro in points.h)
//hx, hy stepsizes
long double update(struct point p, struct point *A[Ny],long double hx, long double hy)
{
	int i;
	int j;
	int row;
	int col;
	long double s;

	long double U;
	long double U_east;
	long double U_west;
	long double U_north;
	long double U_south;

	//Store U value of point and coordinates of point
	row = p.y;
	col = p.x;
	s = p.s;
	
	//printf("Updating A[%d][%d] \n",row,col);	
	
	// Find neighbors of point, only use if in the grid
	// Assigning INFTY ensures the neighbor won't be used in update
	// calculation
	if (row > 0 )
	{	
		U_north = A[row-1][col].U;
		//printf("Assigned north as A[%d][%d].U=%f \n",row-1,col,A[row-1][col].U);
	}
	else
	{
		U_north = INFTY;
	}
	if (row < Ny -1)
	{
		U_south= A[row+1][col].U;
		//printf("Assigned south as A[%d][%d].U=%f \n",row+1,col,A[row+1][col].U);
	}
	else
	{
		U_south = INFTY;
	}

	if (col > 0)
	{
		U_west =A[row][col-1].U;
		//printf("Assigned west as A[%d][%d].U=%f \n",row,col-1,A[row][col-1].U);

	}
	else
	{
		U_west = INFTY;
	}

	if (col < Nx -1)
	{
		U_east = A[row][col+1].U;
		//printf("Assigned east as A[%d][%d].U=%f \n",row,col+1,A[row][col+1].U);
	}
	else
	{
		U_east = INFTY;
	}

	if (U_south == INFTY && U_west == INFTY && U_north == INFTY  && U_east == INFTY)
	{
		printf("Error: All Neighbors have U = infty \n");
	}

	//values for calculating quadratic for 2 point update
	long double a;
	long double b;
	long double c;
	long double roots[2];

	//U value and neighboring update values
	long double U_new;
	long double U_ne;
	long double U_nw;
	long double U_se;
	long double U_sw;

	//////Solve for NE quadrant(See Vladimirsky paper)
	a = hx*hx + hy*hy;
	b = -2*(hy*hy*U_east + hx*hx*U_north);
	c = hy*hy*U_east*U_east + hx*hx*U_north*U_north - hx*hx*hy*hy*s*s;
	solve_quadratic(a,b,c, &roots[0]);

	//2nd root is always larger, don't need to take max of roots
	if (roots[1] >= fmax(U_east,U_north))
	{
		U_ne = roots[1];
	}
	else if (U_east < U_north)
	{
		U_ne = U_east + hx*s;
	}
	else
	{
		U_ne = U_north + hy*s;
	}
	
	/////Solve for NW quadrant
	a = hx*hx + hy*hy;
	b = -2*(hy*hy*U_west + hx*hx*U_north);
	c = hy*hy*U_west*U_west + hx*hx*U_north*U_north - hx*hx*hy*hy*s*s;
	solve_quadratic(a,b,c, &roots[0]);
	if (roots[1] >= fmax(U_west,U_north))
	{
		U_nw = roots[1];
	}
	else if (U_west < U_north)
	{
		U_nw = U_west + hx*s;
	}
	else
	{
		U_nw = U_north + hy*s;
	}
	
	////Solve for SE Quadrant
	a = hx*hx + hy*hy;
	b = -2*(hy*hy*U_east + hx*hx*U_south);
	c = hy*hy*U_east*U_east + hx*hx*U_south*U_south - hx*hx*hy*hy*s*s;
	solve_quadratic(a,b,c, &roots[0]);
	if (roots[1] >= fmax(U_east,U_south))
	{
		U_se = roots[1];	
	}
	else if (U_east < U_south)
	{
		U_se = U_east + hx*s;
	}
	else
	{
		U_se = U_south + hy*s;
	}
	////Solve for SW quadrant
	a = hx*hx + hy*hy;
	b = -2*(hy*hy*U_west + hx*hx*U_south);
	c = hy*hy*U_west*U_west + hx*hx*U_south*U_south - hx*hx*hy*hy*s*s;
	solve_quadratic(a,b,c, &roots[0]);
	if (roots[1] >= fmax(U_west,U_south))
	{
		U_sw = roots[1];	
	}
	else if (U_west < U_south)
	{
		U_sw = U_west + hx*s;
	}
	else
	{
		U_sw = U_south + hy*s;
	}
	
	////find minimum of quadrant U values
	U_new = fmin(U_ne,fmin(U_nw,fmin(U_sw,U_se)));
	//printf("New U: %f \n",U_new);	
	return U_new;
	
	}

