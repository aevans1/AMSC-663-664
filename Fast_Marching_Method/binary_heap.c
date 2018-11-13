#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//Code for implementing a binary heap structure for floats
//Oct. 13th, 2018

//Print heap level by level, k-level heap with N elements
void print_heap(float *heap, int k, int N)
{
	int i;
	int j;
	int q;

	for (i = 0; i < k; i++) 
	{
		printf("Level %d \n",i);
		
		j = 0;
		q = pow(2,i) + j - 1;
		while (q < N)
		{
			printf("heap[%d]= %f \n",q,heap[q]);
			j++;
		}
		printf("\n");

	}
}

//Swaps values of two float variables
void swap_float(float *p, float *q)
{
		float temp;
		temp = *p;
		*p = *q;
		*q = temp;
}

//Add new value to end of heap, re-balances heap
//p is value added to heap
//count is number of last inex of heap
int add_heap(float *heap, float p, int *count)
{
	int i;

	//printf("adding %f \n",p);
	heap[*count] = p;

	if (*count > 0)
	{
		float temp_parent; // Float for checking on negative values of parent indices
		int parent;
		int current;

		//Start at end node, calculate parent node
		current = *count;
		temp_parent = (current-1)/2.0;
	
		//Move new node up in heap until rebalanced
		while (temp_parent >= 0)
		{
			 parent = (int)temp_parent; //round and use as parent index
			
			/*DEBUGGING
			printf("parent after rounding is = %d \n",parent);
			printf("heap[%d]= %f\n",parent,heap[parent]);
			printf("Current heap vals: \n");
				
			for (i = 0; i <= *count; i++) 
			{
				//printf("heap[%d]= %f \n",i,heap[i]);
			}
			*/

			//Swap parent, child if child smaller value	
			//Otherwise, add point to end of heap
			if (heap[current] < heap[parent])	
			{
				/*DEBUGGING*/	
				//printf("In func, heap[%d] = %f \n",parent,heap[parent]);
				//printf("In func, heap[%d] = %f \n",current,heap[*count]);
				//printf("Above value should be %f \n",p);

				swap_float(&heap[current],&heap[parent]);
				
				/*DEBUGGING*/
				//printf("Swapped to...\n");
				//printf("In func, heap[%d] = %f \n",parent,heap[parent]);
				//printf("In func, heap[%d] = %f \n",current,heap[current]);
		
			//Update child, parent values
			temp_parent = (parent -1)/2.0;
			current = parent;

			}
			else
			{
				//Break the loop, Heap is re-balanced
				temp_parent = -1;
			}
		}
	}
	//printf("exited loop for adding %f \n",p);
	//printf("\n");
	
	(*count)++;

	return 0;
}

float pop_heap(float *heap, int *count, int k)
{
	//top of the heap	
	float root;
	root = heap[0];
	
	//printf("popping: %f \n",root);
	//printf("end of heap: %f \n",heap[*count - 1]);

	//Replace top of heap with end, delete end of heap
	swap_float(&heap[0],&heap[*count - 1]);
	heap[*count - 1] = -1000;	
	//printf("new root: %f \n",heap[0]);
	(*count)--;
	
	int lchild;
	int rchild;
	int current;
	int min_child;

	
	current = 0;
	lchild = 1;
	
	//Move new root down in heap until rebalanced	
	while(lchild < *count)
	{
		/*Debugging ////
		//printf("Current: heap[%d] = %f \n",current,heap[current]);
		//printf("Left child: heap[%d] = %f \n",lchild,heap[lchild]);
		*/

		//Assign right child, if right child is in the array
		if (2*current + 1 < *count - 1)
			{	
				rchild = 2*current + 2;
				//printf("Right child: heap[%d] = %f \n \n",rchild,heap[rchild]);
			}
		else
			{
				rchild = lchild;
			}
		
		//Pick min of children
		min_child = rchild;
		if (heap[lchild]<heap[rchild])
		{
			min_child = lchild;
		}	
	
		//Swap current with child if unbalanced
		if (heap[current] > heap[min_child] )
		{
			swap_float(&heap[current],&heap[min_child]);
			//printf("Swapping \n");	
			current = min_child;
			lchild = 2*current + 1;

			/*DEBUGGING, Print Heap
			int q;
			int w;
			int j;
	
			for (q = 0; w < k; w++) 
			{
				printf("Level %d \n",w);
				
				j = 0;
				q = pow(2,w) + j - 1;
				while (q < *count)
				{
					printf("heap[%d]= %f \n",q,heap[q]);
					j++;
				}
				printf("\n");
			}
			*/

		}
	
		//Otherwise, break the loop, heap is re-balanced
		else
		{
			lchild = *count;
		}
	}
	//return popped root value
	return root;
}

int main()
{
	int i;
	int j;
	int k;
	int N;
	k = 3; //number of tree levels
	N = pow(2,k) - 1; //number of nodes in tree, k = 1 is root

	float my_heap[N];
	float my_floats[N];

//Initialize Random float array
	srand( time(NULL) );
	
	for (i = 0; i < N;i++)
	{	
		my_floats[i] = rand();
	}

	int count;
	count = 0; //index of end of heap
	
//Filling heap	
	for (i = 0; i < N; i++)
	{
		add_heap(&my_heap[0],my_floats[i],&count);
	}

// Testing add_heap
	int q;
	int w;
	for (i = 0; i < k; i++) 
	{
		printf("Level %d \n",i);
		for (j = 0; j < pow(2,i) && j < count; j++)
		{
			q = pow(2,i) + j - 1;	
			if (q < count)
				{
					printf("heap[%d]= %f \n",q,my_heap[q]);
				}
		}
		printf("\n");
	}

//Testing pop_heap
	while (count > 0)
	{
	pop_heap(&my_heap[0],&count,k);
	}


//TESTING for heap property
	printf("\n\n Testing for heap property: \n");
	int parent;
	for (i = 0; i < count; i++)
	{
		parent = (i - 1)/2;
		if (my_heap[i] < my_heap[parent])
		{
			printf("error! \n");
			printf("heap[%d] = %f \n",i,my_heap[i]);
			printf("heap[%d] = %f \n",parent,my_heap[parent]);
			printf("\n");
		}
	}
}

