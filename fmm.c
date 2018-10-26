#include <stdio.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>
#include "fmm.h"

//TODO: comments here



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

int main()
{
	int i,j;

	//Define stepsizes, heap for Trial U values, count for iteration number
	long double hx, hy;
	hx = (XMAX - XMIN)/(Nx - 1);
	hy = (YMAX - YMIN)/(Ny - 1);	

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
			A[i][j].x = XMIN + hy*j;
			A[i][j].y = YMIN + hx*i;
			A[i][j].s = 1;
			//A[i][j].s = 1.0/(2.0 + 5.0*A[i][j].x + 20.0*A[i][j].y);
			//printf("s = %f \n",A[i][j].s);
			A[i][j].U = INFTY;
		}
	}

	struct point *heap;
	heap = (struct point*)malloc(Nx*Ny*sizeof(struct point));

	int count;
	count = 0;


	//Define initial boundary Gamma 
	//Using two point sources:
	long double init[2][2]; //global coordinates of initial boundary in R^2
	long double temp;	
	int row,col;

	//set initial two points as (0,0) and (0.8,0)
	int num_initial;
	num_initial = 2;
	init[0][0] = 0.0;
	init[0][1] = 0.0;
	init[1][0] = 0.8;
	init[1][1] = 0.0;
	
	////////////////////////////////
	/*Initialization of algorithm */
	///////////////////////////////
	//Label all neighbors of Known points as Trial, update U values of Neighbors and add their U-values to heap
	
	int new_row, new_col;
	long double h;	

	for (i = 0; i < num_initial; i++)
	{
		//Find A array indices for inital coordinates, label as Known
		temp = round((init[i][0] - XMIN)/hy);
		row = temp;
		temp = round((init[i][1] - YMIN)/hx);
		col = temp;	

		printf("Row: %d, Col: %d \n",row,col);
		A[row][col].label = '2';
		A[row][col].U = 0.0;
	}
	
	for (i = 0; i < num_initial; i++)
	{
		temp = round((init[i][0] - XMIN)/hy);
		row = temp;
		temp = round((init[i][1] - YMIN)/hx);
		col = temp;	

		for (j = 0; j < 4; j++)
		{
			int	neighbor[4][2]= {{row+1, col},{row-1,col},{row,col+1},{row,col-1}};	

			new_row = neighbor[j][0];
			new_col = neighbor[j][1];
	
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
	long double temp_update;

	//continue until heap is empty	
	while (count > 0)
	{
		//find point with lowest U value, label as known
		new_known = pop_heap(&heap[0],&count);
		
		temp = round((new_known.x - XMIN)/hy);
		row = temp;
		temp = round((new_known.y - YMIN)/hx);
		col = temp;	
		
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
	long double tmp,tmp1,tmp2,aux_x,aux_y,max_err = 0,err,s;
	for (i = 0; i < Ny; i++)
	{
		aux_y = YMIN + hx*i;
		for (j = 0; j < Nx; j++)
		{
			fprintf(fid,"%.6e\t",A[i][j].U);
			aux_x = XMIN + hy*j;
			s = 1.0/(2.0 + 5.0*aux_x + 20.0*aux_y);
			tmp1 = (1.0/sqrt(425.0))*acosh(1.0 + 0.5*0.5*s*425.0*((aux_x - 0)*(aux_x-0) + (aux_y - 0)*(aux_y - 0)));
			tmp2 = (1.0/sqrt(425.0))*acosh(1.0 + (1.0/6.0)*0.5*s*425.0*((aux_x - 0.8)*(aux_x-0.8) + (aux_y - 0)*(aux_y - 0)));
			tmp = fmin(tmp1,tmp2);
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



