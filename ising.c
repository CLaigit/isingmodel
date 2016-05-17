#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */

#define  C_a 1.0
#define  MASS 1.0
#define  OMEGA 0.15
#define  EPSILON C_a
#define  NTIMESTEPS 1000
#define  WARMSTEPS 100000
#define  POSSWEEPS 1000000
#define  POSMEASURE 100
#define  DELTA 3
#define  ROW (NTIMESTEPS + 2)
#define  COLUMN (POSSWEEPS/POSMEASURE)
#define  NOUT (1.0 * POSSWEEPS/POSMEASURE)
#define  MIU2 2
#define  LAMBDA 1.0

double kinetic(double x1, double x2){
    return 0.5 * MASS * (x1 - x2) * (x1 - x2) / EPSILON / EPSILON;
}

double potential(double x1, double x2){
    double meanx = 0.5 * (x1 + x2);
    return 0.5 * MASS * OMEGA * OMEGA * meanx * meanx;
}

double totalenergy(double x1, double x2){
    return kinetic(x1, x2) + potential(x1, x2);
}

int main (int argc, char *argv[]){

    static double pos[COLUMN][ROW] = {};
    double shift;
    double deltas;
    double tmppos[ROW] = {};
    static double energy[COLUMN] = {};
    static double exppectPos[COLUMN] = {};
    int warm_memory;
    double updatex;

    srand (time(NULL));
    // Start from to calculate the second column
    for (int i = 1; i < NOUT; i++){
        // Copy the position to an tempory position array
        for(int j = 0; j < ROW; j++){
            tmppos[j] = pos[i-1][j];
        }
        // Check if it is the warmup process or the regular process
        warm_memory = ((i == 1)? WARMSTEPS: POSMEASURE);
        // Run warmup steps or regular steps
        for (int k = 0; k < warm_memory; k++){
            //The first element stores the last position
            // The last element stores the first position
            tmppos[0] = tmppos[ROW-2];
            tmppos[ROW-1] = tmppos[1];
            // Update every point in a path
            for(int t = 1; t <= NTIMESTEPS; t++){
                updatex = tmppos[t] + 2.0 * DELTA * ((double)rand() / (double)RAND_MAX - 0.5);
                deltas = EPSILON * totalenergy(tmppos[t+1], updatex);
                deltas += EPSILON * totalenergy(updatex, tmppos[t-1]);
                deltas -= EPSILON * totalenergy(tmppos[t+1], tmppos[t]);
                deltas -= EPSILON * totalenergy(tmppos[t], tmppos[t-1]);

                if (deltas < 0 || (double)rand() / (double)RAND_MAX <= exp(-deltas)){
                    tmppos[t] = updatex;
                }
            }
            //Update the boundary condition
            tmppos[0] = tmppos[ROW-2];
            tmppos[ROW-1] = tmppos[1];
        }
        // Copy the tempory array to the position matrix
        for(int j = 0; j < ROW; j++){
            pos[i][j] = tmppos[j];
        }
    }
    // Output data
    for(int i = 1; i < COLUMN; i++){
        for(int j = 1; j < ROW-2; j++){
            printf("%f,", pos[i][j]);
        }
        printf("%f\n", pos[i][ROW-2]);
    }
}
