/*
Ising model: Halmitonian H = /sum_ij J(sigma_i)(sigma_j)
We set J = 1 first
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */

#define  LATTICE_LENGTH 200
#define  COLUMN LATTICE_LENGTH
#define  ROW LATTICE_LENGTH
#define  BOLTZMANN_CONST 1
#define  WARMSTEPS 1e3
#define  NSWEEPS 100

// Calculate the energy of the (up, center) (down, center) (left, center) ( right, center)
double energy(int up, int down, int left, int right, int center){
    return -center * (up + down + left + right);
}

int main (int argc, char *argv[]){


    static int lattice[LATTICE_LENGTH][LATTICE_LENGTH] = {};
    double T = 2;
    int col, row;

    T = argc > 1 ? atof(argv[1]) : 2;
    col = argc > 2 ? atoi(argv[2]) : 20;
    row = col;
    // Tempurature
    int new;
    double beta = 1.0 / BOLTZMANN_CONST / T;
    double deltaE = 0.0;
    double tmpE = 0.0, tmpE2 = 0.0, averE = 0.0, averE2 = 0.0;
    double tmpmag = 0.0, tmpmag2 = 0.0, avermag = 0.0, avermag2 = 0.0;
    double siteE = 0.0;
    srand (time(NULL));
    // Initialize every grid point
    for (int i = 0; i < COLUMN; i++){
        for(int j = 0; j < ROW; j++){
            lattice[i][j] = 2 * (rand() % 2) - 1;
            // lattice[i][j] = 1;
        }
    }
    // Warmup process
    for (int nstep = 0; nstep < WARMSTEPS; nstep++){
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
    for (int nstep = 0; nstep < NSWEEPS; nstep++){
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
        tmpE = 0, tmpE2 = 0, tmpmag = 0, tmpmag2 = 0;
        for(int i = 0; i < COLUMN; i++){
            for(int j = 0; j < ROW; j++){
                siteE += energy(lattice[ (i - 1 + ROW) % ROW][j], lattice[(i + 1 + ROW) % ROW][j], lattice[i][(j - 1 + ROW) % ROW], lattice[i][(j + 1 + ROW) % ROW], lattice[i][j])/2;
                tmpE += siteE;
                tmpE2 += siteE * siteE;
                tmpmag += lattice[i][j];
                tmpmag2 += lattice[i][j] * lattice[i][j];
            }
        }
        averE += (1.0 * tmpE / LATTICE_LENGTH / LATTICE_LENGTH / NSWEEPS);
        averE2 += (1.0 * tmpE2 / LATTICE_LENGTH / LATTICE_LENGTH / NSWEEPS);
        avermag += (1.0 * tmpmag / LATTICE_LENGTH / LATTICE_LENGTH / NSWEEPS);
        avermag2 += (1.0 * tmpmag2 / LATTICE_LENGTH / LATTICE_LENGTH / NSWEEPS);
        // Output data every NOUT
        // for(int i = 0; i < COLUMN; i++){
        //     for(int j = 0; j < COLUMN-1; j++){
        //         printf("%d,", lattice[i][j]);
        //     }
        //     printf("%d\n", lattice[i][COLUMN-1]);
        // }
    }
    printf("%f\n", T);
    printf("%f\n", averE);
    printf("%f\n", 1.0*(averE2 - averE * averE) / T / T);
    printf("%f\n", 1.0*(avermag2 - avermag * avermag) / T );
}
