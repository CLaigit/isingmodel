/*
Ising model: Halmitonian H = /sum_ij J(sigma_i)(sigma_j)
We set J = 1 first

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */

#define  LATTICE_LENGTH 20
#define  COLUMN LATTICE_LENGTH
#define  ROW LATTICE_LENGTH
#define  BOLTZMANN_CONST 1
#define  WARMSTEPS 1e6
#define  MEASURESTEPS 1e6
#define  NOUT 1e2
#define  NSWEEPS (MEASURESTEPS/NOUT)

// Calculate the energy of the (up, center) (down, center) (left, center) ( right, center)
int energy(int up, int down, int left, int right, int center){
    double H;
    H = -up * center;
    H -= down * center;
    H -= left * center;
    H -= right * center;
    return H;
}

int main (int argc, char *argv[]){

    static int lattice[COLUMN][ROW] = {};
    double T = 2;

    T = argc > 1 ? atof(argv[1]) : 2;
    // Tempurature
    int new;
    double beta = 1.0 / BOLTZMANN_CONST / T;
    double deltaE;
    srand (time(NULL));
    // Initialize every grid point
    for (int i = 0; i < COLUMN; i++){
        for(int j = 0; j < ROW; j++){
            lattice[i][j] = 2 * (rand() % 2) - 1;
        }
    }
    // Warmup process
    for (int inter = 0; inter < WARMSTEPS; inter++){
        for(int i = 0; i < COLUMN; i++){
            for(int j = 0; j < ROW; j++){
                // flip a spin
                // If the energy becomes small, accept it
                // Else accept w.r.t the probability e^-beta*deltaE
                new = -lattice[i][j];
                deltaE = energy(lattice[ (i - 1 + ROW) % ROW][j], lattice[(i + 1 + ROW) % ROW][j], lattice[i][(j - 1 + ROW) % ROW], lattice[i][(j + 1 + ROW) % ROW], new);
                deltaE -= energy(lattice[ (i - 1 + ROW) % ROW][j], lattice[(i + 1 + ROW) % ROW][j], lattice[i][(j - 1 + ROW) % ROW], lattice[i][(j + 1 + ROW) % ROW], lattice[i][j]);
                if ((double)rand() / (double)RAND_MAX <= exp(- beta * deltaE)){
                    lattice[i][j] = new;
                }
            }
        }
    }

    // Measure steps
    for (int inter = 0; inter < NSWEEPS; inter++){
        // for(int k = 0; k < NOUT; k++){
            for(int i = 0; i < COLUMN; i++){
                for(int j = 0; j < ROW; j++){
                    new = -lattice[i][j];
                    deltaE = energy(lattice[ (i - 1 + ROW) % ROW][j], lattice[(i + 1 + ROW) % ROW][j], lattice[i][(j - 1 + ROW) % ROW], lattice[i][(j + 1 + ROW) % ROW], new);
                    deltaE -= energy(lattice[ (i - 1 + ROW) % ROW][j], lattice[(i + 1 + ROW) % ROW][j], lattice[i][(j - 1 + ROW) % ROW], lattice[i][(j + 1 + ROW) % ROW], lattice[i][j]);
                    if (deltaE < 0 || (double)rand() / (double)RAND_MAX <= exp(- beta * deltaE)){
                        lattice[i][j] = new;
                    }
                }
            }
        // }
        // Output data every NOUT
        for(int i = 0; i < COLUMN; i++){
            for(int j = 0; j < COLUMN-1; j++){
                printf("%d,", lattice[i][j]);
            }
            printf("%d\n", lattice[i][COLUMN-1]);
        }
    }
}
