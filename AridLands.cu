/*
ROHIT GUPTA & JOHAN VAN DE KOPPEL
Arid Pattern formation
June 2010
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// include CUDA
#include <cuda_runtime.h>

// include parameter kernel
#include "AridLands.h"

////////////////////////////////////////////////////////////////////////////////
// Allocates a matrix with random float entries
////////////////////////////////////////////////////////////////////////////////

void randomInit(float* data, int x_siz, int y_siz, int type)
{
	int i,j;
	for(i=0;i<y_siz;i++)
	{
		for(j=0;j<x_siz;j++)
		{
			//for every element find the correct initial
			//value using the conditions below
			if(i==0||i==y_siz-1||j==0||j==x_siz-1)
			    data[i*x_siz+j]=0.0f; // This value for the boundaries
			else
			{
				if(type==Plants)
					
				{
                    // A randomized initiation here
					if((rand() / (float)RAND_MAX)<0.05f)
							data[i*x_siz+j] = 100.0f;
						else
							data[i*x_siz+j] = 0.0f;
				}
				else if(type==Surface_Water)
					data[i*x_siz+j]=(float)(R/(alpha*W0));
				else if(type==Soil_Water)
					data[i*x_siz+j]=(float)(R/rw/4);
			}
		}
	}			

} // End randomInit

////////////////////////////////////////////////////////////////////////////////
// Laplacation operator definition, to calculate diffusive fluxes
////////////////////////////////////////////////////////////////////////////////

__device__ float
LaplacianXY(float* pop, int row, int column)
{
	float retval;
	int current, left, right, top, bottom;	
	float dx = dX;
	float dy = dY;
	
	current=row * WidthGrid + column;	
	left=row * WidthGrid + column-1;
	right=row * WidthGrid + column+1;
	top=(row-1) * WidthGrid + column;
	bottom=(row+1) * WidthGrid + column;

	retval = ( (( pop[current] - pop[left] )/dx ) 
		      -(( pop[right] - pop[current] )/dx )) / dx + 
		     ( (( pop[current] - pop[top] )/dy  ) 
			  -(( pop[bottom] - pop[current] )/dy ) ) / dy;

	return retval;
}

////////////////////////////////////////////////////////////////////////////////
// Gradient operator definition, to calculate advective fluxes
////////////////////////////////////////////////////////////////////////////////


__device__ float
GradientY(float* pop, int row, int column)
{
	float retval;
	int current, top;	
	float dy =dY;
	
	current=row * WidthGrid + column;	
	top=(row-1) * WidthGrid + column;
	
	retval =  (( pop[current] - pop[top] )/dy ); 

	return retval;
}

////////////////////////////////////////////////////////////////////////////////
// Simulation kernel
////////////////////////////////////////////////////////////////////////////////

__global__ void 
Aridlandskernel (float* popP, float* popW, float* popO)
{

	//run for U X V times. For every U times completed store in the array storeA and storeM

	int current;

	float d2popPdxy2, d2popWdxy2, d2popOdxy2;
	float drP,drW, drO;
	
	int row		=	blockIdx.y*Block_Size_Y+threadIdx.y;
	int column	=	blockIdx.x*Block_Size_X+threadIdx.x;
	
	current=row * WidthGrid + column;
	
	if(row > 0 && row < HeightGrid-1 && column > 0 && column < WidthGrid-1)
	 {		
		//Now calculating terms for the O Matrix	
		d2popOdxy2 = -DifO * LaplacianXY(popO, row, column) - AdvO * GradientY(popO, row, column);
		drO = (R-alpha*(popP[current]+k2*W0)/(popP[current]+k2)*popO[current]);
		 
		//Now calculating terms for the W Matrix
		d2popWdxy2 = -DifW * LaplacianXY(popW, row, column);
		drW = (alpha*(popP[current]+k2*W0)/(popP[current]+k2)*popO[current] 
	          - gmax*popW[current]/(popW[current]+k1)*popP[current]-rw*popW[current]);

		//Now calculating terms for the P Matrix
		d2popPdxy2 = -DifP * LaplacianXY(popP, row, column);
		drP = (c*gmax*popW[current]/(popW[current]+k1)*popP[current] - d*popP[current]);

		__syncthreads();

		popO[current]=popO[current]+(drO+d2popOdxy2)*dT;		
		popW[current]=popW[current]+(drW+d2popWdxy2)*dT;	
		popP[current]=popP[current]+(drP+d2popPdxy2)*dT;
	
	 }

	__syncthreads();

	// HANDLE Boundaries
	if(row==0)
			{
				popW[row * WidthGrid + column]=popW[(HeightGrid-2) * WidthGrid+column];
				popO[row * WidthGrid + column]=popO[(HeightGrid-2) * WidthGrid+column];
				popP[row * WidthGrid + column]=popP[(HeightGrid-2) * WidthGrid+column];
			}
	else if(row==HeightGrid-1)			
			{
				popW[row * WidthGrid + column]=popW[1*WidthGrid+column];
				popO[row * WidthGrid + column]=popO[1*WidthGrid+column];
				popP[row * WidthGrid + column]=popP[1*WidthGrid+column];
			}	
	else if(column==0)			
			{
				popW[row * WidthGrid + column]=popW[row * WidthGrid + WidthGrid-2];
				popO[row * WidthGrid + column]=popO[row * WidthGrid + WidthGrid-2];
				popP[row * WidthGrid + column]=popP[row * WidthGrid + WidthGrid-2];
			}	
	else if(column==WidthGrid-1)			
			{
				popW[row * WidthGrid + column]=popW[row * WidthGrid + 1];
				popO[row * WidthGrid + column]=popO[row * WidthGrid + 1];
				popP[row * WidthGrid + column]=popP[row * WidthGrid + 1];
			}	
	
} // End Aridlandskernel

////////////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char** argv)
{
	int store_count;       // The number of times a frame was stored
	double time_elapsed;   // The amount of time that has passed

	int NumStored;
	int store_i,i;
	unsigned int size_grid = WidthGrid * HeightGrid;
	unsigned int mem_size_grid = sizeof(float) * size_grid;
	unsigned int size_storegrid = WidthGrid * HeightGrid * MAX_STORE;
	unsigned int mem_size_storegrid = sizeof(float) * size_storegrid;

    double RecordTimes[MAX_STORE];
    
	float* h_store_popP;
	float* h_store_popO;
	float* h_store_popW;

	float* h_popP;
	float* h_popO;
	float* h_popW;

	float* d_popP;
	float* d_popO;
	float* d_popW;
    
	FILE *fp;
    
	int height_matrix=HeightGrid;
	int width_matrix=WidthGrid;
    
    /*--------------------INITIALIZATIONS ON HOST-------------------*/
	time_elapsed=Time;
	store_count=0;

	// set seed for rand()
	srand((unsigned)time( NULL ));

	//allocate host memory for matrices popP, popO, and popW
	h_popP = (float*) malloc(mem_size_grid);
	h_popO = (float*) malloc(mem_size_grid);
	h_popW = (float*) malloc(mem_size_grid);	

	//allocate host memory for matrices store_popP, store_popO, and store_popW
	h_store_popP=(float*) malloc(mem_size_storegrid);
	h_store_popO=(float*) malloc(mem_size_storegrid);
	h_store_popW=(float*) malloc(mem_size_storegrid);

    /*---------------------INITIALIZING THE ARRAYS----------------------------*/
	randomInit(h_popP, WidthGrid, HeightGrid, Plants);
	randomInit(h_popO, WidthGrid, HeightGrid, Surface_Water);
	randomInit(h_popW, WidthGrid, HeightGrid, Soil_Water);

    /*---------------------INITIALIZING ON GPU--------------------------------*/
	
        // allocate device memory    
	cudaMalloc((void**) &d_popP, mem_size_grid);
	cudaMalloc((void**) &d_popO, mem_size_grid);
	cudaMalloc((void**) &d_popW, mem_size_grid);

        //copy host memory to device
	cudaMemcpy(d_popP, h_popP, mem_size_grid, cudaMemcpyHostToDevice);
	cudaMemcpy(d_popO, h_popO, mem_size_grid, cudaMemcpyHostToDevice);
	cudaMemcpy(d_popW, h_popW, mem_size_grid, cudaMemcpyHostToDevice);   
 
    /*---------------------SETUP EXECUTION PARAMETERS-------------------------*/	
	dim3 threads;      // Setting up the GPU setting, thread block size
	dim3 grid;         // Setting up the GPU setting, grid structure

	threads.x= Block_Size_X;
	threads.y= Block_Size_Y;
	grid.x=DIVIDE_INTO(WidthGrid,threads.x);
	grid.y=DIVIDE_INTO(HeightGrid,threads.y);

    //using namespace std;
    clock_t begin = clock();

    // Calculate the times at which the simulation is stored
    for(i=0;i<=NumFrames;i++) 
		{ RecordTimes[i]=(double)i/(double)NumFrames*(double)EndTime; }  

    /*----- Printing info to the screen --------------------------------*/
	system("clear");
        printf("\n");
	printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n");
	printf(" * Arid Land Patters                                     * \n");		
	printf(" * CUDA implementation : Rohit Gupta, 2009               * \n");
	printf(" * Following a model by Rietkerk et al 2002              * \n");
	printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n");
	
	int device;
    cudaGetDevice(&device);

    struct cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);
    
    printf(" Using device: %s \n\n", props.name);

	printf(" Current grid dimensions: %d x %d cells\n\n", WidthGrid, HeightGrid);
    
    /*----- The simulation loop ----------------------------------------*/
	while((float)time_elapsed<=(float)EndTime)
	   {
		// execute the kernel
  		Aridlandskernel<<< grid, threads >>>(d_popP, d_popW, d_popO);

		// Storing the data if a particular number of timesteps have passed	
		if((float)time_elapsed>=(float)RecordTimes[store_count])
		    {		
			cudaMemcpy((void *)h_popP, (void *)d_popP, mem_size_grid,
					   cudaMemcpyDeviceToHost);
			cudaMemcpy((void *)h_popO, (void *)d_popO, mem_size_grid,
					   cudaMemcpyDeviceToHost);	
			cudaMemcpy((void *)h_popW, (void *)d_popW, mem_size_grid,
					   cudaMemcpyDeviceToHost);
			
			//Store values at this frame.
			memcpy(h_store_popP+(store_count*size_grid),h_popP,mem_size_grid);
			memcpy(h_store_popO+(store_count*size_grid),h_popO,mem_size_grid);
			memcpy(h_store_popW+(store_count*size_grid),h_popW,mem_size_grid);

			fprintf(stderr, "\r Current timestep: %1.0f of %1.0d, Storepoint %1d of %1d",
                                        time_elapsed,EndTime, store_count, NumFrames);
	
                        #if SaveEachPlot == on

			//	fp=fopen("CurrentFrame.dat","wb");	

			//	// Storing parameters
			//	fwrite(&width_matrix,sizeof(int),1,fp);
			//	fwrite(&height_matrix,sizeof(int),1,fp);

		        //      fwrite(&h_popP,sizeof(float),size_grid,fp);
                	//	fwrite(&h_popO,sizeof(float),size_grid,fp);
                	//	fwrite(&h_popW,sizeof(float),size_grid,fp);
	
			//	fclose(fp);
                        #endif	

			store_count=store_count+1;		
				
		    } // if on writing one frame ends
	
		time_elapsed=time_elapsed+(double)dT;
	
	   } //while on time ends
	
	/*---------------------Report on time spending----------------------------*/
	// sdkStopTimer(&timer);
	clock_t end = clock();
	//printf("\n\n %1.5f to %1.5f \n", (float)begin, (float)end);
    double elapsed_secs = (double)(end - begin)/CLOCKS_PER_SEC;
	printf("\n\n Processing time: %4.1f (s) \n", elapsed_secs);
	//sdkDeleteTimer(&timer);

	/*---------------------Write to file now----------------------------------*/
	fp=fopen("AridLands.dat","wb");	
	NumStored=store_count;	

	// Storing parameters
	fwrite(&width_matrix,sizeof(int),1,fp);
	fwrite(&height_matrix,sizeof(int),1,fp);
	fwrite(&NumStored,sizeof(int),1,fp);

    for(i=0;i<=NumStored;i++)
	   { fwrite(&RecordTimes[i],sizeof(double),1,fp); }	
	
	for(store_i=0;store_i<store_count;store_i++)
	   {
		fwrite(&h_store_popP[store_i*size_grid],sizeof(float),size_grid,fp);
		fwrite(&h_store_popO[store_i*size_grid],sizeof(float),size_grid,fp);
		fwrite(&h_store_popW[store_i*size_grid],sizeof(float),size_grid,fp);

		printf("\r Saving simulation results: %2.0f%%", (float)store_i/(float)store_count*100.0);
	   }
	
	printf("\r Saving simulation results: 100%%\n\n");

	fclose(fp);

	/*---------------------Clean up memory------------------------------------*/
	free(h_popP);
	free(h_popO);
	free(h_popW);

	free(h_store_popP);
	free(h_store_popO);
	free(h_store_popW);	
    
	cudaFree(d_popP);
	cudaFree(d_popO);
	cudaFree(d_popW);

	cudaThreadExit();
	cudaDeviceReset();
	
	#if defined(__APPLE__) && defined(__MACH__)
    system("say All ready");
    #endif

}
