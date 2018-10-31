/*Each point in domain A has associated struct 'point':
 -label is for FMM labels
 -U is the current U-value at point
 -x,y are coordinates of point
 -Point labels:
	0 - Far
	1 - Trial
	2 - Known
*/
struct point
{
	char label;
	double U;	

	int row;
	int col;
	double s;

};
typedef struct point point;

/*struct for containing (x,y) coordinates in R^2, a 2D vector*/
struct vect
{
	double x;
	double y;
};
typedef struct vect vect;

////Physical width and length of domain

//One point source:
#define XMIN -3.5
#define XMAX 3.5
#define YMIN -0.7
#define YMAX 0.7

//Two point sources:
//#define XMIN 0.0
//#define XMAX 1.0
//#define YMIN 0.0
//#define YMAX 1.0

////number of steps in x,y directions
#define Nx 17
#define Ny 17

////Max possible value for U, default for Far points
#define INFTY 1e6

//Declaring all point_heap functions
void swap_point(struct point *p, struct point *q);
void print_heap (struct point *heap, int k, int N);
void add_heap(struct point *heap, struct point p, int *count);
struct point pop_heap(struct point *heap, int *count);

//Declaring all helper functions
double update(struct point p, struct point *A[Ny], double hx, double hy);
int in_mesh(int row, int col);
void get_coord(int row, int col, double hx, double hy, vect *v);
void get_meshindex(int *row, int *col, double hx, double hy, vect v);
void get_neighbors(int neighbors[4][2],int row, int col);

//Declare quadratic function
void solve_quadratic(double a, double b, double c, double *roots);

