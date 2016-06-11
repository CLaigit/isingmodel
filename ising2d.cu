/*
 *  Ising model: Halmitonian H = /sum_ij J(sigma_i)(sigma_j)
 */

/*
 *   1. Calculate the energy in the program
 *   2. Calculate the heat capacity in the program
 *   3. Add more inputs to adjust the length of lattice
 *   4. A matlab code to plot data.
 *       data format example:
 *                    position.x  position.y   spin(-1, 1)
 *       Iteattion 1:    1           4               -1
 *                       *           *                *
 *                       *           *                *
 *       Iteattion 2:    4           3                1
 *                       *           *                *
 *                       *           *                *
 *       Iteattion N:    35          76               1
 *                       *           *                *
 *                       *           *                *
 *   5. Compare the numerical value with the analytic value
 *   6. Move to 3D
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */
#include <curand.h>
#include <curand_kernel.h>

/*
 * LATTICE_LENGTH  - length of the lattice
 * LATTICE_LENGTH  - number of element is one lattice
 * BOLTZMANN_CONST - bolzmann constant. Set to 1.
 */

#define LATTICE_LENGTH 1024
#define LATTICE_2 (LATTICE_LENGTH * LATTICE_LENGTH)
#define BOLTZMANN_CONST 1
#define N LATTICE_LENGTH

#define WARM_STEP 1e3
#define MEAS_STEP 1e3
#define WARP 1e1
#define NUM_THREAD_X 32
#define NUM_THREAD_Y 32
#define TEMPERATURE 4.0

__device__ int energy(int up, int down, int left, int right, int center);
__global__ void update(int *lattice, double beta, double *E_d, double *M_d, double *E2_d, double *M2_d, int tag, curandState * global_state);
__global__ void printstate(int *lattice);
__global__ void init_rand(curandState * global_state, unsigned long seed);


/* Setup random seed to each kernel */
__global__ void init_rand(curandState * global_state, unsigned long seed){
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	curand_init(seed, idx + idy * N, 0, &global_state[idx + idy * N]);
	__syncthreads();
}

/*
 *   update is the function to update a point
 *   1. flip a point (1 -> -1 or -1 -> 1)
 *   2. compare the energy before flip a point and after flip a point
 *   3. if the energy with flipped point is small, accept
 *   4. if the energy is larger, generate a random number pro_rand (0,1),
 *      if pro_rand < e^(-beta * delatE), aceept. else reject.
 */
__global__ void update(int* lattice, double beta, double *E_d, double *M_d, double *E2_d, double *M2_d, int tag, curandState * global_state){
	// Calculate the global index
	// Calculate the global index for the up, down, left, right index.
		
	// declare parameters
	int itx, ity, idx, idy, index;
	int flip, up, down, left, right, center;
	double pro_rand, deltaE, E;

	// local index
	itx = threadIdx.x;
	ity = threadIdx.y;

	// global index
	idx = blockIdx.x * blockDim.x + itx;
	idy = blockIdx.y * blockDim.y + ity;
	index = idx * N + idy;
		
	// load data into shared memory
	__shared__ int lat[32 + 2][32 + 2];
	__syncthreads();

	lat[itx+1][ity+1] = lattice[index];

	if(idx == 0){
		lat[itx][ity + 1] = lattice[index + (N - 1) * N];
	}else if(itx == 0){
		lat[itx][ity + 1] = lattice[index - N];
	}

	if(idx == N - 1){
		lat[itx + 2][ity + 1] = lattice[index - (N - 1) * N];
	}else if(itx == NUM_THREAD_X - 1){
		lat[itx + 2][ity + 1] = lattice[index + N -1];
	}

	if(idy == 0){
		lat[itx + 1][ity] = lattice[index + N - 1];
	}else if(ity == 0){
		lat[itx + 1][ity] = lattice[index - 1];
	}

	if(idy == N - 1){
		lat[itx + 1][ity + 2] = lattice[index - (N - 1)];
	}else if(ity == NUM_THREAD_X - 1){
		lat[itx + 1][ity + 2] = lattice[index + 1];
	}
		
	curandState local_state = global_state[idx * N + idy];
	pro_rand = curand_uniform(&local_state);
	global_state[idx * N + idy] = local_state;
	__syncthreads();

	// for even sites
	if((idx + idy) % 2 == 0){
   	up     = lat[itx][ity + 1];
   	down   = lat[itx + 2][ity + 1];
   	left   = lat[itx + 1][ity];
   	right  = lat[itx + 1][ity + 2];
   	center = lat[itx + 1][ity + 1];

		// Flip the center element
		flip = -center;

		// Calculate the difference between these two state
		E      = energy(up, down, left, right, center);
		deltaE = -2.0 * E;

		// If deltaE < 0 or pro_rand <= e^(-beta * deltaE), accept new value
		if (deltaE < 0 || pro_rand <= exp(- 1.0 * beta * (deltaE * 1.0))){
      lat[itx + 1][ity + 1] *= -1;
     }
	}

	// wait for even site completion
	__syncthreads();

	// for odd sites
	if((idx + idy) % 2 == 1){
		up     = lat[itx][ity + 1];
    down   = lat[itx + 2][ity + 1];
    left   = lat[itx + 1][ity];
    right  = lat[itx + 1][ity + 2];
    center = lat[itx + 1][ity + 1];
	
		// Flip the center element
		flip = -center;

		// Calculate the difference between these two state
    E      = energy(up, down, left, right, center);
		deltaE = -2.0 * E;

		// If deltaE < 0 or pro_rand <= e^(-beta * deltaE), accept new value
		if (deltaE < 0 || pro_rand <= exp(- 1.0 * beta * (deltaE * 1.0))){
			lat[itx + 1][ity + 1] *= -1;
		}
	}

	// wait for odd site completion
	__syncthreads();
			

	// store data back
	lattice[index] = lat[itx + 1][ity + 1];

	if(tag == 1){
		E_d[index] += E;
		M_d[index] += lat[itx+1][ity+1];
		E2_d[index] += E * E;
		M2_d[index] += lat[itx+1][ity+1] * lat[itx+1][ity+1];
	}
	__syncthreads();

}

/*
 *   printstate is the function to print the whole matrix.
 *   Since it prints in parallel, we also print the global
 *   index of the matrx.
 *   it prints (x, y, (1 or -1)).
 */
__global__ void printstate(int* lattice) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if (idx < N && idy < N){
		printf("%d, %d, %d\n", idx, idy, lattice[idx + idy * N]);
	}
	__syncthreads();
}

/*
 *   energy is the function used to calculate the energy between
 *   (center, up), (center, down), (center, left), (center, right)
 */
__device__ int energy(int up, int down, int left, int right, int center){
	double H;
	H = - up * center - down * center - left * center - right * center;
	return H;
}

/*
 *   Commandline inputs option
 *   1. Tempurature (T)
 *
 */
int main (int argc, char *argv[]){

	int *lattice;
	int *d_lattice;

	double *E;
	double *E_d;

	double *E2;
	double *E2_d;

	double *M;
	double *M_d;

	double *M2;
	double *M2_d;

	double T = TEMPERATURE;
 	int warmsteps = WARM_STEP;
	int nout = MEAS_STEP;
	int warp = WARP;

	int numthreadx = NUM_THREAD_X;
	int numthready = NUM_THREAD_Y;
	int numblocksX = LATTICE_LENGTH / numthreadx;
	int numblocksY = LATTICE_LENGTH / numthready;

	// First input: Tempurature. Usually between (1, 6),
	// Critical Tempurature is around 2.2
	T = argc > 1 ? atof(argv[1]) : T;
	warmsteps = argc > 2 ? atof(argv[2]) : warmsteps;
	nout = argc > 3 ? atof(argv[3]) : nout;
	warp = argc > 4 ? atof(argv[4]) : warp;

	// Define the size of lattice and energy
	const size_t bytes_lattice = LATTICE_2 * sizeof(int);
	const size_t bytes_E = LATTICE_2 * sizeof(double);
	const size_t bytes_M = LATTICE_2 * sizeof(double);

	// Allocate memory for lattice. It is a lattice^2 long array.
	// The value can only be 1 or -1.
	lattice = (int*)malloc(LATTICE_2 * sizeof(int));

	E = (double*)malloc(LATTICE_2 * sizeof(double));
	M = (double*)malloc(LATTICE_2 * sizeof(double));

	E2 = (double*)malloc(LATTICE_2 * sizeof(double));
	M2 = (double*)malloc(LATTICE_2 * sizeof(double));
		

	// initialize lattice by rand(-1, 1)
	for(int i = 0; i < LATTICE_2; i++){
		lattice[i] = 2 * (rand() % 2) - 1;
		E[i] = 0.0;
		M[i] = 0.0;
		E2[i] = 0.0;
		M2[i] = 0.0;
   }

	// Set dimensions of block and grid
	dim3 grid(numblocksX, numblocksY, 1);
	dim3 thread(numthreadx, numthready,1);

	// set up random for each kernel 
	curandState *global_state;
	cudaMalloc(&global_state, LATTICE_2 * sizeof(curandState));
	init_rand<<< grid, thread >>> (global_state, unsigned(time(NULL)));

	// beta is a parameter in the probability
	double beta = 1.0 / (BOLTZMANN_CONST * 1.0) / T;

	// Allocate memoery in device and copy from host to device
	cudaMalloc((void **)&d_lattice, bytes_lattice);
	cudaMalloc((void **)&E_d, bytes_E);
	cudaMalloc((void **)&M_d, bytes_M);

	cudaMalloc((void **)&E2_d, bytes_E);
	cudaMalloc((void **)&M2_d, bytes_M);

	cudaMemcpy(d_lattice, lattice, bytes_lattice, cudaMemcpyHostToDevice);
	cudaMemcpy(E_d, E, bytes_E, cudaMemcpyHostToDevice);
	cudaMemcpy(M_d, M, bytes_M, cudaMemcpyHostToDevice);

	cudaMemcpy(E2_d, E2, bytes_E, cudaMemcpyHostToDevice);
	cudaMemcpy(M2_d, M2, bytes_M, cudaMemcpyHostToDevice);

	// To change the buffer size of printf; otherwise it cannot print all data
	cudaDeviceSetLimit(cudaLimitPrintfFifoSize, N * N * sizeof(int));

//	printf("Testing for T = %2f, beta = %2f...\n", T, beta);
	
	// Warmup process
//	printf("Starting Warming Steps... \n");
	int cnt = 0;

	for (int iter = 0; iter < warmsteps; iter++){
//		printf("\r [ %f% ] ", (100.0 * cnt++) / warmsteps);
			
		update<<<grid, thread>>>(d_lattice, beta, E_d, M_d, E2_d, M2_d, 0, global_state);
		cudaDeviceSynchronize();
	}
//	printf("\n");

	// Measure process
//	printf("Starting Measurement Steps... \n");
	cnt = 0;
	int cnt2 = 0;

	for (int nstep = 0; nstep < nout; nstep++){
//		printf("\r [ %f% ] ", (100.0 * cnt++) / nout);

		if(nstep % warp == 0){
			cnt2++;
     	update<<<grid, thread>>>(d_lattice, beta, E_d, M_d, E2_d, M2_d, 1, global_state);
		}else{
     	update<<<grid, thread>>>(d_lattice, beta, E_d, M_d, E2_d, M2_d, 0, global_state);
		}
		cudaDeviceSynchronize();

	}
//	printf("\n");
		
	double energy = 0.0;
	double magnetization = 0.0;

	double energy2 = 0.0;
	double magnetization2 = 0.0;

	cudaMemcpy(lattice, d_lattice, bytes_E, cudaMemcpyDeviceToHost);

	cudaMemcpy(E, E_d, bytes_E, cudaMemcpyDeviceToHost);
	cudaMemcpy(M, M_d, bytes_M, cudaMemcpyDeviceToHost);

	cudaMemcpy(E2, E2_d, bytes_E, cudaMemcpyDeviceToHost);
	cudaMemcpy(M2, M2_d, bytes_M, cudaMemcpyDeviceToHost);
		
	for(int i = 0; i < LATTICE_2; i++){
		energy += E[i];
		magnetization += M[i];

		energy2 += E2[i];
		magnetization2 += M2[i];
	}
	
	double avg_E = energy / cnt2 / (LATTICE_2 * 1.0) / 2.0;
	double avg_M = magnetization / cnt2 / (LATTICE_2 * 1.0);
	avg_M = avg_M < 0 ? -avg_M : avg_M;

	double avg_E2 = energy2 / cnt2 / (LATTICE_2 * 1.0) / 4.0;
	double avg_M2 = magnetization2 / cnt2 / (LATTICE_2 * 1.0);

	double heat_cap = 1.0 * (avg_E2 - avg_E * avg_E) / T / T;
	double mag_sus  = 1.0 * (avg_M2 - avg_M * avg_M) / T; 

//	printf("Average energy: %5f \n", avg_E);
//	printf("Average magnetization: %5f \n", avg_M);
	printf("%5f %5f %5f %5f %5f\n", T, avg_E, avg_M, heat_cap, mag_sus);

	free(lattice);
	free(E);
	free(M);
	free(E2);
	free(M2);
	cudaFree(d_lattice);
	cudaFree(E_d);
	cudaFree(M_d);
	cudaFree(E2_d);
	cudaFree(M2_d);
		
	return 0;
}
