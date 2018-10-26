/*Each point in domain A has associated struct 'point':
 -label is for FMM labels
 -U is the current U-value at point
 -x,y are coordinates of point
 -Point labels:
	0 - Far
	1 - Trial
	2 - Accepted
*/
struct point
{
	char label;
	long double U;	

	long double x;
	long double y;
	long double s;

};

////Physical width and length of domain
#define XMIN 0.0
#define XMAX 1.0
#define YMIN 0.0
#define YMAX 1.0

////number of steps in x,y directions
#define Nx 257
#define Ny 257

////Max possible value for U, default for Far points
#define INFTY 1e6

//Declaring all point_heap functions
void swap_point(struct point *p, struct point *q);
void print_heap (struct point *heap, int k, int N);
void add_heap(struct point *heap, struct point p, int *count);
struct point pop_heap(struct point *heap, int *count);

//Declaring update function
long double update(struct point p, struct point **A, long double hx, long double hy);

//Declare quadratic function
void solve_quadratic(long double a, long double b, long double c, long double *roots);

