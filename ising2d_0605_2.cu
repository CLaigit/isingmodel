/*
 *  Ising model: Halmitonian H = /sum_ij J(sigma_i)(sigma_j)
 */

/*
 * TODO:
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
//#include <stdint.h>

/*
 * LATTICE_LENGTH is the length of the lattice
 * LATTICE_LENGTH is the number of element is one lattice
 * BOLTZMANN_CONST is bolzmann constant. It is set to 1.
 */

#define LATTICE_LENGTH 256
#define LATTICE_2 (LATTICE_LENGTH * LATTICE_LENGTH)
#define BOLTZMANN_CONST 1
#define N LATTICE_LENGTH

#define WARM_STEP 1e4
#define MEAS_STEP 1e4
#define WARP 1e2
#define NUM_THREAD_X 16
#define NUM_THREAD_Y 16

__device__ int energy(int up, int down, int left, int right, int center);
__global__ void update(int *lattice, double beta, double *rand_prob, double *ttl_E, int tag);
__global__ void printstate(int *lattice);


/*
 *   update is the function to update a point
 *   1. flip a point (1 -> -1 or -1 -> 1)
 *   2. compare the energy before flip a point and after flip a point
 *   3. if the energy with flipped point is small, accept
 *   4. if the energy is larger, generate a random number pro_rand (0,1),
 *      if pro_rand < e^(-beta * delatE), aceept. else reject.
 */
__global__ void update(int* lattice, double beta, double *rand_prob, double *ttl_E, int tag){
    // Calculate the global index
    // Calculate the global index for the up, down, left, right index.
		
		// declare parameters
		int itx, ity, idx, idy, index;
    int flip, up, down, left, right, center;
    double pro_rand;
    double deltaE;
		double E;

		// index in the block
		itx = threadIdx.x;
		ity = threadIdx.y;

		// global index
    idx = blockIdx.x * blockDim.x + itx;
    idy = blockIdx.y * blockDim.y + ity;
		index = idx + idy * N;
		
		// load data into shared memory
		__shared__ int lat[16 + 2][16 + 2];
		__syncthreads();

		lat[itx+1][ity+1] = lattice[index];

		if(idx == 0){
			lat[itx][ity + 1] = lattice[index + N - 1];
		}else if(itx == 0){
			lat[itx][ity + 1] = lattice[index - 1];
		}

		if(idx == N - 1){
			lat[itx + 2][ity + 1] = lattice[index - N + 1];
		}else if(itx == 16 - 1){
			lat[itx + 2][ity + 1] = lattice[index + 1];
		}

		if(idy == 0){
			lat[itx + 1][ity] = lattice[index + (N-1) * N];
		}else if(ity == 0){
			lat[itx + 1][ity] = lattice[index - N];
		}

		if(idy == N - 1){
			lat[itx + 1][ity + 2] = lattice[index - (N-1) * N];
		}else if(ity == 16-1){
			lat[itx + 1][ity + 2] = lattice[index + N];
		}
		
		pro_rand = rand_prob[index];	
		__syncthreads();

		// for even sites
		if(index % 2 == 0){
    	up     = lat[itx][ity + 1];
    	down   = lat[itx + 2][ity + 1];
    	left   = lat[itx + 1][ity];
    	right  = lat[itx + 1][ity + 2];
    	center = lat[itx + 1][ity + 1];

      // Flip the center element
      flip = -center;

      // Calculate the difference between these two state
			deltaE = - flip * (up + down + left + right);
			E      = - center * (up + down + left + right);
			deltaE -= E;
//      deltaE = energy(up, down, left, right, flip);
//      deltaE -= energy(up, down, left, right, center);

      // If deltaE < 0 or pro_rand <= e^(-beta * deltaE), accept new value
      if (deltaE < 0 || pro_rand <= exp(- beta * deltaE)){
      	lat[itx + 1][ity + 1] = flip;
      }
		}

		// wait for even site completion
		__syncthreads();

		// for odd sites
		if(index % 2 == 1){

    	up     = lat[itx][ity + 1];
    	down   = lat[itx + 2][ity + 1];
    	left   = lat[itx + 1][ity];
    	right  = lat[itx + 1][ity + 2];
    	center = lat[itx + 1][ity + 1];

      // Flip the center element
      flip = -center;

      // Calculate the difference between these two state
			deltaE = - flip * (up + down + left + right);
			E      = - center * (up + down + left + right);
			deltaE -= E;
//      deltaE = energy(up, down, left, right, flip);
//      deltaE -= energy(up, down, left, right, center);

      // If deltaE < 0 or pro_rand <= e^(-beta * deltaE), accept new value
      if (deltaE < 0 || pro_rand <= exp(- beta * deltaE)){
      	lat[itx + 1][ity + 1] = flip;
      }
		}

		// wait for odd site completion
		__syncthreads();
			

		// store data back
		lattice[index] = lat[itx + 1][ity + 1];
//		printf("%d, %d, %d\n", idx, idy, lattice[index]);		
		if(tag == 1){
			ttl_E[index] += E;
//printf("----->>>>> cur_energy %5f\n", E);
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
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

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
    H = -up * center;
    H -= down * center;
    H -= left * center;
    H -= right * center;
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
		double *rand_prob;
		double *rand_prob_d;
		double *E;
		double *ttl_E;

    double T = 2.0;
    int warmsteps = WARM_STEP;
    int nout = MEAS_STEP;
    int warp = WARP;

    int numthreadx = NUM_THREAD_X;
    int numthready = NUM_THREAD_Y;
    int numblocksX = LATTICE_LENGTH / numthreadx;
    int numblocksY = LATTICE_LENGTH / numthready;

    // First input: Tempurature. Usually between (1, 6),
    // Critical Tempurature is around 2.2
    T = argc > 1 ? atof(argv[1]) : 2;

    // Define the size of lattice and energy
    const size_t bytes_lattice = LATTICE_2 * sizeof(int);
		const size_t bytes_rand = LATTICE_2 * sizeof(double);
		const size_t bytes_E = LATTICE_2 * sizeof(double);

    // Allocate memory for lattice. It is a lattice^2 long array.
    // The value can only be 1 or -1.
    lattice = (int*)malloc(LATTICE_2 * sizeof(int));
		rand_prob = (double*)malloc(LATTICE_2 * sizeof(double));
		E = (double*)malloc(LATTICE_2 * sizeof(double));
    // energy = (int*)malloc(sizeof(int));
		

		srand(time(NULL));
    // initialize lattice by rand(-1, 1)
    for(int i = 0; i < LATTICE_2; i++){
			lattice[i] = 2 * (rand() % 2) - 1;
			E[i] = 0.0;
    }

    // Set dimensions of block and grid
    dim3 grid(numblocksX, numblocksY, 1);
    dim3 thread(numthreadx, numthready,1);

    // beta is a parameter in the probability
    double beta = 1.0 / (BOLTZMANN_CONST * 1.0) / T;

    // Allocate memoery in device and copy from host to device
    cudaMalloc((void **)&d_lattice, bytes_lattice);
		cudaMalloc((void **)&rand_prob_d, bytes_rand);
		cudaMalloc((void **)&ttl_E, bytes_E);
    // cudaMalloc((void **)&d_energy, bytes_energy);

    cudaMemcpy(d_lattice, lattice, bytes_lattice, cudaMemcpyHostToDevice);
    cudaMemcpy(ttl_E, E, bytes_E, cudaMemcpyHostToDevice);
    // cudaMemcpy(d_energy, energy, bytes_energy, cudaMemcpyHostToDevice);
//    cudaMemcpy(ttl_E, (double)0.0, sizeof(double), cudaMemcpyHostToDevice);

    // To change the buffer size of printf; otherwise it cannot print all data
    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, N * N * sizeof(int));

    // Warmup process
printf("Starting Warming Steps... \n");
int cnt = 0;
    for (int iter = 0; iter < warmsteps; iter++){
printf("\r [ %f% ] ", (100.0 * cnt++) / warmsteps);
			
    	for(int i = 0; i < LATTICE_2; i++){
				rand_prob[i] = ((double)rand() / (double)RAND_MAX);
            //lattice[i] = 1;
    	}
			cudaMemcpy(rand_prob_d, rand_prob, bytes_rand, cudaMemcpyHostToDevice);
      cudaDeviceSynchronize();

      update<<<grid, thread>>>(d_lattice, beta, rand_prob_d, ttl_E, 0);
      cudaDeviceSynchronize();
    }
printf("\n");

    // Measure process
printf("Starting Measurement Steps... \n");
cnt = 0;
    for (int nstep = 0; nstep < nout; nstep++){
printf("\r [ %f% ] ", (100.0 * cnt++) / nout);

    	for(int i = 0; i < LATTICE_2; i++){
				rand_prob[i] = ((double)rand() / (double)RAND_MAX);
            //lattice[i] = 1;
    	}
			cudaMemcpy(rand_prob_d, rand_prob, bytes_rand, cudaMemcpyHostToDevice);
      cudaDeviceSynchronize();

//      update<<<grid, thread>>>(d_lattice, beta, rand_prob_d, ttl_E, 0);
			if(nout % warp == 0){
      	update<<<grid, thread>>>(d_lattice, beta, rand_prob_d, ttl_E, 1);
			}else{
      	update<<<grid, thread>>>(d_lattice, beta, rand_prob_d, ttl_E, 0);
			}
      cudaDeviceSynchronize();

//      if(nstep % warp == 0){
//	      printstate<<<grid, thread>>>(d_lattice);
//      	cudaDeviceSynchronize();
//			}
    }
printf("\n");
		
		double energy = 0;
		cudaMemcpy(E, ttl_E, bytes_E, cudaMemcpyDeviceToHost);
		cudaMemcpy(lattice, d_lattice, bytes_E, cudaMemcpyDeviceToHost);
		
		for(int i = 0; i < LATTICE_2; i++){
			energy += E[i];
		}
	
		double avg_E = energy / (nout * 1.0 / warp) / (LATTICE_2 * 1.0);
		printf("acerage energy: %5f \n", avg_E);

    free(lattice);
    cudaFree(d_lattice);
		free(rand_prob);
		cudaFree(rand_prob_d);
		cudaFree(ttl_E);
}
