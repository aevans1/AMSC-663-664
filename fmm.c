#include <stdio.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>
#include "fmm.h"


////////////////////////////////////////////
//Fast Marching Method for solving the eikonal equation ||grad U|| = s(x)
//To change parameters, see header file fmm.h
//To compile, type 'make' in command line
///////////////////////////////////////////
int main()
{
	int i,j;
	double hx,hy;
	vect v;
	//int Nx,Ny;
	//Nx = 257;	
	//Ny = 257;
	//Define stepsizes, heap for Trial U values, count for iteration number
	hx = (XMAX - XMIN)/(Nx - 1.0);
	hy = (YMAX - YMIN)/(Ny - 1.0);	

	printf("hx = %f \n",hx);	
	//Initialize Domain
	
	//point *A;	
	//A = (point*)malloc(Nx*Ny*sizeof(point));	
	
	point  *A[Ny];
	for (i = 0; i < Ny; i++)
	{
		A[i] = (point *)malloc(Ny*sizeof(point));
	}

	double x,y;
	for (i = 0; i < Ny; i++)
	{
		for (j = 0; j < Nx; j++)
		{
			A[i][j].label = '0'; //label all as 'Far'
			A[i][j].row = i;
			A[i][j].col = j;

			//Two point sources, specific speed function
			get_coord(i,j,hx,hy,&v);		
			A[i][j].s = 1.0/(2.0 + 5.0*v.x + 20.0*v.y);
			//A[i][j].s = 1; //speed function identically 1, 1 point source
			A[i][j].U = INFTY;
		}
	}
	point *heap;
	heap = (point*)malloc(Nx*Ny*sizeof(point));

	int count;
	count = 0;

	//Define initial boundary
	//
	//Using two point sources:

	//Define point source
	//Mark center as 'Known', set U to 0

	//make array of initial coords for boundary
	int num_initial;
	
	/*One point source*/
	//num_initial = 1;
	//vect init[num_initial];
	//init[0].x = 0.0;
	//init[0].y = 0.0;

	/*Two point sources(see Cameron's note)*/
	num_initial = 2;
	vect init[num_initial];
	init[0].x = 0.0;
	init[0].y = 0.0;
	init[1].x = 0.8;
	init[1].y = 0.0;

	printf("init[1] = \n");
	printf("%f \n",init[1].x);
	printf("%f \n",init[1].y);

	////////////////////////////////
	/*Initialization of algorithm */
	///////////////////////////////
	
	//Label all neighbors of Known points as Trial, update U values of Neighbors and add their U-values to heap
	int row, col, new_row, new_col;
	double h;
	int neighbor[4][2];

	for(i = 0; i < num_initial; i++)
	{
		get_meshindex(&row, &col, hx, hy, init[i]);
		
		A[row][col].label = '2';
		A[row][col].U = 0.0;
	}

	for(i = 0; i < num_initial; i++)
	{
		get_meshindex(&row, &col, hx, hy, init[i]);
		//int	neighbor[4][2]= {{row+1, col},{row-1,col},{row,col+1},{row,col-1}};	
		get_neighbors(neighbor,row,col);

		for (j = 0; j < 4; j++)
		{
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
	point new_known;
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
		//int neighbor[4][2] = {{row+1, col},{row-1,col},{row,col+1},{row,col-1}};	
		get_neighbors(neighbor,row,col);
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
	FILE *gid;
	gid = fopen("trueU.txt","w");
	double tmp,aux_x,aux_y,max_err = 0,err,tmp1,tmp2,s;
	for (i = 0; i < Ny; i++)
	{
		aux_y = YMIN + hy*i;
		for (j = 0; j < Nx; j++)
		{
 			//printf("%0.2f\n ",A[i][j].U);
			//printf("%0.2f \n",*A(i*Nx + j).U);
			fprintf(fid,"%.6e\t",A[i][j].U);
			aux_x = XMIN + hx*j;
			//tmp = sqrt(aux_x*aux_x + aux_y*aux_y);
			s = 1.0/(2.0 + 5.0*aux_x + 20.0*aux_y);
			tmp1 = (1.0/sqrt(425.0))*acosh(1.0 + 0.5*0.5*s*425.0*((aux_x - 0)*(aux_x-0) + (aux_y - 0)*(aux_y - 0)));
			tmp2 = (1.0/sqrt(425.0))*acosh(1.0 + (1.0/6.0)*0.5*s*425.0*((aux_x - 0.8)*(aux_x-0.8) + (aux_y - 0)*(aux_y - 0)));
			tmp = fmin(tmp1,tmp2);
			fprintf(gid,"%.6e\t",tmp);
			err = A[i][j].U - tmp;
			if( err > max_err ) max_err = err;
		}
// 		printf("\n");
		fprintf(fid,"\n");
		fprintf(gid,"\n");
	}
	fclose(fid);
	fclose(gid);
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



