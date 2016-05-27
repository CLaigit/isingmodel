/*
Ising model: Halmitonian H = /sum_ij J(sigma_i)(sigma_j)
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

/*
* LATTICE_LENGTH is the length of the lattice
* LATTICE_LENGTH is the number of element is one lattice
* BOLTZMANN_CONST is bolzmann constant. It is set to 1.
*/

#define  LATTICE_LENGTH 256
#define  LATTICE_2 (LATTICE_LENGTH * LATTICE_LENGTH)
#define  BOLTZMANN_CONST 1
#define  N LATTICE_LENGTH

__device__ int energy(int up, int down, int left, int right, int center);
__global__ void update(int *lattice, unsigned int offset);
__global__ void printstate(int *lattice);


/*
*   update is the function to update a point
*   1. flip a point (1 -> -1 or -1 -> 1)
*   2. compare the energy before flip a point and after flip a point
*   3. if the energy with flipped point is small, accept
*   4. if the energy is larger, generate a random number pro_rand (0,1),
*      if pro_rand < e^(-beta * delatE), aceept. else reject.
*/
__global__ void update(int* lattice, const unsigned int offset, double beta){
    // Calculate the global index
    // Calculate the global index for the up, down, left, right index.
    const unsigned int idx = blockIdx.x * blockDim.y + threadIdx.x + offset;
    const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y + offset;
    const unsigned int idx_l = (idx - 1 + N) % N;
    const unsigned int idx_r = (idx + 1 + N) % N;
    const unsigned int idy_u = (idy - 1 + N) % N;
    const unsigned int idy_d = (idy + 1 + N) % N;
    int flip, up, down, left, right, center;
    double pro_rand;
    double deltaE;

    // To generate random number in cuda
    curandState_t state;
    curand_init(idx, idx + 1, 0, &state);

    if (idx < N && idy < N && idx_l < N && idx_r < N && idy_u < N && idy_d < N){
        // generate a random number between (0,1) uniformly
        pro_rand = curand_uniform(&state);

        up = lattice[idx + idy_u * N];
        down = lattice[idx + idy_d * N];
        left = lattice[idx_l + idy * N];
        right = lattice[idx_r + idy * N];
        center = lattice[idx + idy * N];

        // Flip the center element
        flip = -center;
        // Calculate the difference between these two state
        deltaE = energy(up, down, left, right, flip);
        deltaE -= energy(up, down, left, right, center);

        // If deltaE < 0 or pro_rand <= e^(-beta * deltaE), accept new value
        if (pro_rand <= exp(- beta * deltaE)){
            lattice[idx + idy * N ] = flip;
        }
    }
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
    double T = 2;
    int warmsteps = 1;
    int nout;
    nout = 2;
    int warp = 1000;

    int numthreadx = 16;
    int numthready = 16;
    int numblocksX = LATTICE_LENGTH / numthreadx;
    int numblocksY = LATTICE_LENGTH / numthready;

    // First input: Tempurature. Usually between (1, 6),
    // Critical Tempurature is around 2.2
    T = argc > 1 ? atof(argv[1]) : 2;

    // Define the size of lattice and energy
    const size_t bytes_lattice = LATTICE_2 * sizeof(int);
    const size_t bytes_energy = sizeof(int);

    // Allocate memory for lattice. It is a lattice^2 long array.
    // The value can only be 1 or -1.
    lattice = (int*)malloc(LATTICE_2 * sizeof(int));
    // energy = (int*)malloc(sizeof(int));

    // initialize lattice by rand(-1, 1)
    for(int i = 0; i < LATTICE_2; i++){
            lattice[i] = 2 * (rand() % 2) - 1;
            //lattice[i] = 1;
    }

    // Set dimensions of block and grid
    dim3 grid(numblocksX, numblocksY, 1);
    dim3 thread(numthreadx, numthready,1);

    // beta is a parameter in the probability
    double beta = 1.0 / BOLTZMANN_CONST / T;

    // Allocate memoery in device and copy from host to device
    cudaMalloc((void **)&d_lattice, bytes_lattice);
    // cudaMalloc((void **)&d_energy, bytes_energy);
    cudaMemcpy(d_lattice, lattice, bytes_lattice, cudaMemcpyHostToDevice);
    // cudaMemcpy(d_energy, energy, bytes_energy, cudaMemcpyHostToDevice);

    // To change the buffer size of printf; otherwise it cannot print all data
    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, N * N * sizeof(int) * N);

    // Warmup process
    for (int iter = 0; iter < warmsteps; iter++){
        update<<<grid, thread>>>(d_lattice, 0, beta);
        update<<<grid, thread>>>(d_lattice, 1, beta);
        cudaDeviceSynchronize();
        if(iter % warp == 0)
            fprintf(stderr,"Warmup Iteration: %d\n", iter);
    }

    // Measure process
    for (int nstep = 0; nstep < nout; nstep++){
        update<<<grid, thread>>>(d_lattice, 0, beta);
        update<<<grid, thread>>>(d_lattice, 1, beta);
        cudaDeviceSynchronize();
        printstate<<<grid, thread>>>(d_lattice);
        if(nstep % warp == 0)
            fprintf(stderr,"Measure Iteration: %d\n", nstep);
    }

    free(lattice);
    cudaFree(d_lattice);
}