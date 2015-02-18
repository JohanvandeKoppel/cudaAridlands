/* 
Rohit Gupta
parts of code borrowed/inspired/copied from cambridge course on cuda
http://www.many-core.group.cam.ac.uk/course/
*/

/* Try to keep grid sise a factor of 30 = number of SMs on the processor */
#ifndef _ARIDLANDS_H_
#define _ARIDLANDS_H_

// Compiler directives
#define on              1
#define off             0

#define LoadCurrentPlot off
#define SaveEachPlot    off

// Thread block size
#define Block_Size_X	32               // 16	 
#define Block_Size_Y	32               // 16

// Number of blox
/* I define the Block_Number_ensions of the matrix as product of two numbers
Makes it easier to keep them a multiple of something (16, 32) when using CUDA*/
#define Block_Number_X	16               // 32
#define Block_Number_Y	16               // 32

// Matrix Block_Number_ensions
// (chosen as multiples of the thread block size for simplicity)
#define WidthGrid  (Block_Size_X * Block_Number_X)			// Matrix A width
#define HeightGrid (Block_Size_Y * Block_Number_Y)			// Matrix A height

// DIVIDE_INTO(x/y) for integers, used to determine # of blocks/warps etc.
#define DIVIDE_INTO(x,y) (((x) + (y) - 1)/(y))

// Definition of spatial parameters
#define dX		0.5              // 0.5   The size of each grid cell in X direction
#define dY		0.5              // 0.5   The size of each grid cell in Y direction 

// Process parameters            Original value   Explanation and Units
#define	R		1.05             // 1.0   Rainfall(mm.d-1), original values 0.75, 1.00, and 1.25 
#define alpha	0.2              // 0.2   Prop of surface water (d-1)
#define W0		0.2              // 0.2   Bare soil infiltration rate (-)
#define rw   	0.2	         	 // 0.2   Soil water loss rate due to seepage and evaporation (d-1)
#define c		10	        	 // 10    Plant uptake constant (g.mm-1.m-2)
#define gmax	0.05             // 0.05  Plant growth constant (mm.g-1.m2.d-1)   
#define d		0.25             // 0.25  Plant senescence rate (d-1)
#define k1		5	        	 // 5     Half saturation constant for plant uptake and growth (mm)
#define k2		5	        	 // 5     Half saturation constant for water infiltration (g.m-2)

#define DifP		0.1	         // 0.1   The diffusion constant of P
#define	DifW		0.1	         // 0.1   The diffusion constant of W
#define	DifO		100	         // 100   The diffusion constant of O
#define AdvO    	0            // 10    The advection constant of O

#define Time		0            // 1     Start time of the simulation
#define dT		    0.0005       // 0.0005 Time step
#define NumFrames	300          // 300   Number of times the data is stored
#define MAX_STORE	(NumFrames+2)// Determines the size of the storage array
#define EndTime		1000         // 1000  End time of the simulation

// Name definitions
#define Plants	        101
#define Surface_Water	102
#define	Soil_Water		103

#define HORIZONTAL		201
#define VERTICAL		202

#endif // _ARIDLANDS_H_

