#include <stdio.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>
#include "fmm.h"

/* Checks if array entry A[row][col] is in bounds of mesh Nx x Ny */
int in_mesh(int row, int col)
{
	int in_mesh = 0;	
	if (row >= 0 && row < Ny && col >=  0 && col < Nx)
	{
		in_mesh = 1;
	}
	return in_mesh;
}

struct vect
{
	double x;
	double y;
};

void get_coord(int row, int col, double hx, double hy, struct vect *v, struct point *A[Ny])
{
	v->x = XMIN + col*hx;
	v->y = YMIN + row*hy;
}

void get_meshindex(int *row, int *col, double hx, double hy, struct vect v)
{
	*row = (v.y - YMIN)/hy;
	*col = (v.x - XMIN)/hx;
}

int main()
{
	int i,j;
	
	//Define stepsizes, heap for Trial U values, count for iteration number
	double hx, hy;
	hx = (XMAX - XMIN)/(Nx - 1);
	hy = (YMAX - YMIN)/(Ny - 1);	

	//TESTING:
	struct vect v;

	//Initialize Domain
	struct point  *A[Ny];
	for (i = 0; i < Ny; i++)
	{
		A[i] = (struct point *)malloc(Ny*sizeof(struct point));
	}

	for (i = 0; i < Ny; i++)
	{
		for (j = 0; j < Nx; j++)
		{
			A[i][j].label = '0'; //label all as 'Far'
			A[i][j].s = 1; //speed function identically 1 currently
			A[i][j].row = i;
			A[i][j].col = j;
			//get_coord(A[i][j].row,A[i][j].col,hx,hy,&v,A);
			//printf("v.x: %f v.y: %f \n",v.x,v.y);
			A[i][j].U = INFTY;
		}
	}
	struct point *heap;
	heap = (struct point*)malloc(Nx*Ny*sizeof(struct point));

	int count;
	count = 0;


	//Define initial boundary
	//
	//Using two point sources:

	//Define point source
	//Mark center as 'Known', set U to 0
	
	
//TESTING	
	//int istart = Nx/2, jstart = Ny/2;
	//A[istart][jstart].label = '2';
	//A[istart][jstart].U = 0;

	//make array of initial coords for boundary
	int num_initial;
	num_initial = 1;
	struct vect init[num_initial];
	init[0].x = 0.0;
	init[0].y = 0.0;

	////////////////////////////////
	/*Initialization of algorithm */
	///////////////////////////////
	
	//Label all neighbors of Known points as Trial, update U values of Neighbors and add their U-values to heap
	//int neighbor[4][2] = {{istart+1, jstart},{istart-1,jstart},{istart,jstart+1},{istart,jstart-1}};	
	int row, col, new_row, new_col;
	double h;
	//int neighbor[4][2];

	for(i = 0; i < num_initial; i++)
	{
		get_meshindex(&row, &col, hx, hy, init[i]);
		
		A[row][col].label = '2';
		A[row][col].U = 0.0;
	}

	for(i = 0; i < num_initial; i++)
	{
		get_meshindex(&row, &col, hx, hy, init[i]);
	
		for (j = 0; j < 4; j++)
		{
			int	neighbor[4][2]= {{row+1, col},{row-1,col},{row,col+1},{row,col-1}};	
			
			new_row = neighbor[i][0];
			new_col = neighbor[i][1];

			//use hy if row difference, hx if column difference
			h = fabs( (new_row - row)*hy + (new_col - col)*hx );

			//Change neighbor of Known point to Trial Point, update value and add to
			//heap
			if (in_mesh(new_row,new_col))
			{
				A[new_row][new_col].label = '1';	
				A[new_row][new_col].U = A[row][col].U + h*A[new_row][new_col].s;		
				add_heap(&heap[0],A[new_row][new_col],&count);
			}
		}
	}	

	//////////////////
	/*Main Loop     */
	//////////////////
	//Continue labelling Known points, Update Trial points, Searching for lowest U-values until mesh is done
	struct point new_known;
	//int row;
	//int col;
	double temp_update;

	//continue until heap is empty	
	while (count > 0)
	{
		//find point with lowest U value, label as known
		new_known = pop_heap(&heap[0],&count);
		row = new_known.row;
		col = new_known.col;
		
		//printf("new_known: A[%d][%d] = %f \n",row,col,new_known.U);	
		A[row][col].label = '2';

		//Find all not Known neighbors of 'New Known', label as trial, and
		//update
		int neighbor[4][2] = {{row+1, col},{row-1,col},{row,col+1},{row,col-1}};	
		int new_row, new_col;
		for (i = 0; i < 4; i++)
		{
			new_row = neighbor[i][0];
			new_col = neighbor[i][1];

			//Check if neighbor is in the mesh, then update	
			if (in_mesh(new_row,new_col))
			{
				temp_update = update(A[new_row][new_col],A,hx,hy);

				//only update if it decreases the U-value
				if (temp_update < A[new_row][new_col].U)
				{
					A[new_row][new_col].U = temp_update;
				}
		
				//If a Far point, label as Trial and add to heap
				if (A[new_row][new_col].label == '0')
				{
					A[new_row][new_col].label = '1';
					add_heap(&heap[0],A[new_row][new_col],&count);
				}
			}
		}
	}
	// end main loop
	
	////Print the domain, check U values
	
	FILE *fid;
	fid = fopen("U.txt","w");
	double tmp,aux_x,aux_y,max_err = 0,err;
	for (i = 0; i < Ny; i++)
	{
		aux_y = YMIN + hy*i;
		for (j = 0; j < Nx; j++)
		{
// 			printf("%0.2f ",A[i][j].U);
			fprintf(fid,"%.6e\t",A[i][j].U);
			aux_x = XMIN + hx*j;
			tmp = sqrt(aux_x*aux_x + aux_y*aux_y);
			err = A[i][j].U - tmp;
			if( err > max_err ) max_err = err;
		}
// 		printf("\n");
		fprintf(fid,"\n");
	}
	fclose(fid);
	printf("Nx = %i, Ny = %i, MaxErr = %.4e\n",Nx,Ny,max_err);
	printf("%i\t%i\t%.4e\n",Nx,Ny,max_err);

	
	/*Free up memory*/
	for (i = 0; i < Ny; i++)
	{
		free(A[i]);
	}
	
	free(heap);
}
//End program



