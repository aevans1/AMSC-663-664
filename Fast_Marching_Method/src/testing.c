#include <stdio.h>

int myfun(int *arr,int m, int n)
{
	int i;
	int j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("my_array[%d][%d]=%d \n",i,j, *((arr +(i*n)) + j));	
		}
	}
	return 0;
}

void myotherfun(int* x)
{
	int val = x[1];
	x = &val;
}
int main()
{

	int h[2];
	int g[2];
	g[0] = 1;
	h[0] = g[0];
	printf("a[0] = %d \n",h[0]);
	g[0] = 5;
	printf("a[0] = %d \n",h[0]);




	float y;
	int l;
	y = l/2.0;
	printf("y=%f",y);
	
	int a[2];
	a[0] = 1;
	a[1] = 9;
	myotherfun(&a[0]);
	printf("a[%d] = %d \n \n",0,a[0]);
	printf("a[%d] = %d \n \n",1,a[1]);

	int my_array[3][3];
	int i;
	int j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			my_array[i][j] = 2*i+j;	
		}
	}

// Printing out array in different ways
// Print printing out a[...]
//
//
//	for (i = 0; i < 9; i++)
//	{
//		printf("my_array[%d]=%d \n",i,my_array[i]);
//	}
////Print out a[...]a[...]
//	for (i = 0; i < 3; i++)
//	{
//		for (j = 0; j < 3; j++)
//		{
//			printf("my_array[%d][%d]=%d \n",i,j,my_array[i][j]);	
//		}
//	}
//
	myfun( (int *)my_array,3,3);
}

//End program
