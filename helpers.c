#include <stdio.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>
#include "fmm.h"

//////////////////////////////////////////////////////////
/*Contains various helper functions for main file fmm.c */
///////////////////////////////////////////////////////////

/* Checks if array entry A[row][col] is in bounds of Ny x Nx mesh */
//Note: currently uses global variables Nx,Ny, defined in header fmm.h


void get_neighbors(int neighbors[4][2],int row, int col)
{
	neighbors[0][0] = row + 1;
	neighbors[0][1] = col;
	
	neighbors[1][0] = row - 1;
	neighbors[1][1] = col;
	
	neighbors[2][0] = row;
	neighbors[2][1] = col + 1;
	
	neighbors[3][0] = row;
	neighbors[3][1] = col - 1;
}

int in_mesh(int row, int col,int Nx, int Ny)
{
	int in_mesh = 0;	
	if (row >= 0 && row < Ny && col >=  0 && col < Nx)
	{
		in_mesh = 1;
	}
	return in_mesh;
}
/*given a mesh point with row,col number, returns x,y coords in input vector v*/
////hx, hy are space steps in x,y direction
void get_coord(int row, int col, double hx, double hy, vect *v)
{
	v->x = XMIN + col*hx;
	v->y = YMIN + row*hy;
}

/*given a vector with x,y coords, returns mesh row, col vals in 'row','col'
  pointers*/
void get_meshindex(int *row, int *col, double hx, double hy, vect v)
{
	*row = (v.y - YMIN)/hy;
	*col = (v.x - XMIN)/hx;
}

/*Find neighbors of point p in domain A Compute neighbors and eikonal update */
//input point p to update,
//array of points A, stepsizes hx,hy
//Note: Ny is a global variable in points.h
double update(point p, point *A, double hx, double hy, int Nx, int Ny)
{
	int i;
	int j;
	int row;
	int col;
	double s;

	double U;
	double U_east;
	double U_west;
	double U_north;
	double U_south;

	//Store U value of point and coordinates of point
	row = p.row;
	col = p.col;
	s = p.s;
	
	//printf("Updating A[%d][%d] \n",row,col);	
	
	// Find neighbors of point, only use if in the grid
	// Assigning INFTY ensures the neighbor won't be used in update
	// calculation
	if (row > 0 )
	{	
		U_north = A[(row-1)*Nx + col].U;
		//printf("Assigned north as A[%d][%d].U=%f \n",row-1,col,A[row-1][col].U);
	}
	else
	{
		U_north = INFTY;
	}
	if (row < Ny -1)
	{
		U_south = A[(row+1)*Nx + col].U;
		//printf("Assigned south as A[%d][%d].U=%f \n",row+1,col,A[row+1][col].U);
	}
	else
	{
		U_south = INFTY;
	}

	if (col > 0)
	{
		U_west =A[row*Nx + (col-1)].U;
		//printf("Assigned west as A[%d][%d].U=%f \n",row,col-1,A[row][col-1].U);
	}
	else
	{
		U_west = INFTY;
	}

	if (col < Nx -1)
	{
		U_east = A[row*Nx + (col+1)].U;
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
	double a;
	double b;
	double c;
	double roots[2];

	//U value and neighboring update values
	double U_new;
	double U_ne;
	double U_nw;
	double U_se;
	double U_sw;

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

