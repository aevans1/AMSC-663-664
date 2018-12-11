# Fast Marching Method
This set of code is designed for implementing the Fast Marching Method in 2D. Currently the code will run FMM for any Eikonal equation with point source intial boundary and any user-specified speed function. 

Quick Description:


Editable files for various FMM implementations:
<br/>
Headerfile: 'fmm.h', problem domain defined here, all methods declared here <br/>
Main file: 'fmm.c', the FMM algorithm is run here, speed function and boundary are defined here. <br/>
<br/>
Helper files for above:<br/>
'binary_heap.c': implemntation of binary heap data structure, uses struct 'point' defined in headerfile<br/>
'quadratic.c': simple implementation of quadratic formula<br/>
'helpers.c': contains various helper functions, most importantly, 'update', which defines the upwind finite-difference step needed in 'fmm.c'

## Getting Started

To run everything, run the file 'fmm.c' is the main file which runs FMM. At the current state of thecode, you will need to comment and un-comment to implement examples. For your own speed functions, you'll have to edit the .s function in line 37, and for different point sources you'll edit lines 58-70.

Example 1: One Point Source

For this example, the domain is [-3.5,3.5] x [-0.7, 0.7] and the speed function is s(x) = 1 over the whole domain. The point source is by default the origin, but any point source can be input. 

To set the domain, go to fmm.h and edit the following:

////Physical width and length of domain <br/>
 31 <br/>
 32 //One point source:<br/>
 33 //#define XMIN -3.5<br/>
 34 //#define XMAX 3.5<br/>
 35 //#define YMIN -0.7<br/>
 36 //#define YMAX 0.7<br/>
 37 <br/> 
though this is the default, the code should run with any domain picked.

Uncomment the following lines:

/*One point source example*/ <br/>
 64     //source_x = 0.0; <br/>
 65     //source_y = 0.0; <br/>
 66     //num_initial = 1; <br/>
 67     //vect init[num_initial]; <br/>
 68     //init[0].x = source_x; <br/>
 69     //init[0].y = source_y; <br/>
 
and input (source_x,source_y) for a given point source. 

After, comment out any other lines defining point sources, particularly
  
71     /*Two point source example*/<br/>
72     num_initial = 2;<br/>
73     vect init[num_initial];<br/>
74     init[0].x = 0.0;<br/>
75     init[0].y = 0.0;<br/>
76     init[1].x = 0.8;<br/>
77     init[1].y = 0.0;<br/>

To verify your solution, uncomment

/*For one point source, identical speed*/<br/>   
//tmp = sqrt((aux_x - source_x)*(aux_x-source_x) + (aux_y-source_y)*(aux_y-source    _y);

The tmp variable above gives the true solution U(x) = ||x - source||^2.

Also, be sure to comment out

187             /*For 2 point sources */<br/>
188             s = 1.0/(2.0 + 5.0*aux_x + 20.0*aux_y);<br/>
189             tmp1 = (1.0/sqrt(425.0))*acosh(1.0 + 0.5*0.5*s*425.0*((aux_x - 0)*(aux_x-0) + (au    x_y - 0)*(aux_y - 0)));<br/>
190             tmp2 = (1.0/sqrt(425.0))*acosh(1.0 + (1.0/6.0)*0.5*s*425.0*((aux_x - 0.8)*(aux_x-    0.8) + (aux_y - 0)
(aux_y - 0)));<br/>
191             tmp = fmin(tmp1,tmp2);<br/>
192             //printf("%0.2f ",tmp);<br/>

Before running fmm.c, check the main() lines in fmm.c to ensure you have the desired NX,Ny steps.

### Prerequisites
gcc

For any compiler debugging check the 'makefile' for the current gcc commands listed.

