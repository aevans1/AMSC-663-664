#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "fmm.h"

//Code for implementing a binary heap structure for struct 'points'
//Oct. 13th, 2018

void swap_point(point *p, point *q)
{
        point temp;
        temp = *p;
        *p = *q;
        *q = temp;
}


//Print heap level by level, k-level heap with N elements
void print_heap (point *heap, int k, int N)
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
			printf("heap[%d]= %f \n",q,heap[q].U);
			j++;
		}
		printf("\n");

	}
}
//Add new value to end of heap, re-balances heap
//point is new value being added to heap
//count is final index of heap after adding p
void add_heap(point *heap, point p, int *count)
{
	int i;
	//printf("adding %f \n",p);
	heap[*count] = p;

	//only check for parents if there are parents in the heap!
	if (*count > 0)
	{
		long double temp_parent; // long double for checking on negative values of parent indices
		int parent;
		int current;

		//Start at end node, calculate its parent
		current = *count;
		temp_parent = (current-1)/2.0;
	
		//Move new node up in heap until rebalanced
		while (temp_parent >= 0)
		{
			 parent = (int)temp_parent; //round and use as parent index
			
			/*DEBUGGING/
			printf("parent after rounding is = %d \n",parent);
			printf("heap[%d]= %f\n",parent,heap[parent].U);
			printf("Current heap vals: \n");
				
			for (i = 0; i <= *count; i++) 
			{
				//printf("heap[%d]= %f \n",i,heap[i].U);
			}
			*/////////

			//Swap parent, child if child smaller value	
			//Otherwise, add point to end of heap
			if (heap[current].U < heap[parent].U)	
			{
				/*DEBUGGING*/	
				//printf("In func, heap[%d] = %f \n",parent,heap[parent].U);
				//printf("In func, heap[%d] = %f \n",current,heap[*count].U);
				//printf("Above value should be %f \n",p);
				////////////
				
				swap_point(&heap[current],&heap[parent]);
				
				/*DEBUGGING*/
				//printf("Swapped to...\n");
				//printf("In func, heap[%d] = %f \n",parent,heap[parent].U);
				//printf("In func, heap[%d] = %f \n",current,heap[current].U);
				////////////
				
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
}

//Remove top of heap, re-balances heap
//point is new value being added to heap
//count is final index of heap before rebalancing heap
point pop_heap(point *heap, int *count)
{
	//top of the heap	
	point root;
	root = heap[0];
	
	//printf("popping: %f \n",root);
	//printf("end of heap: %f \n",heap[*count - 1].U);

	//Replace top of heap with end, delete end of heap
	swap_point(&heap[0],&heap[*count - 1]);
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
		/*DEBUGGING
		//printf("Current: heap[%d] = %f \n",current,heap[current].U);
		//printf("Left child: heap[%d] = %f \n",lchild,heap[lchild].U);
		*//////////

		//Assign right child, if right child is in the array
		if (2*current + 1 < *count - 1)
			{	
				rchild = 2*current + 2;
				//printf("Right child: heap[%d] = %f \n \n",rchild,heap[rchild].U);
			}
		else
			{
				rchild = lchild;
			}
		
		//Pick min of children
		min_child = rchild;
		if (heap[lchild].U<heap[rchild].U)
		{
			min_child = lchild;
		}	
	
		//Swap current with child if unbalanced
		if (heap[current].U > heap[min_child].U )
		{
			swap_point(&heap[current],&heap[min_child]);
			//printf("Swapping \n");	
			current = min_child;
			lchild = 2*current + 1;
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

/*Main loop for testing
//int main()
{
	//int i;
	//int j;
	//int k;
	//int N;
	//k = 3; //number of tree levels
	//N = pow(2,k) - 1; //number of nodes in tree, k = 1 is root

	//long double my_heap[N];
	//long double my_long doubles[N];

//In//itialize Random long double array
	//srand( time(NULL) );
	//
	//for (i = 0; i < N;i++)
	//{	
	//	my_long doubles[i] = rand();
	//}

	//int count;
	//count = 0; //index of end of heap
	//
//Fi//lling heap	
	//for (i = 0; i < N; i++)
	//{
	//	add_heap(&my_heap[0],my_long doubles[i],&count);
	//}

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

//Te//sting pop_heap
	//while (count > 0)
	//{
	//pop_heap(&my_heap[0],&count,k);
	//}


//TESTING for heap property
	printf("\n\n Testing for heap property: \n");
	int parent;
	for (i = 0; i < count; i++)
	{
		parent = (i - 1)/2;
		if (my_heap[i].U < my_heap[parent].U)
		{
			printf("error! \n");
			printf("heap[%d].U = %f \n",i,my_heap[i].U);
			printf("heap[%d].U = %f \n",parent,my_heap[parent].U);
			printf("\n");
		}
	}
	return 0;
}
*/

