#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>       /* time */
#define LATTICE_LENGTH 100

int main (int argc, char *argv[]){

    static int lattice[LATTICE_LENGTH] = {};
    srand (time(NULL));

    for (int i = 0; i < LATTICE_LENGTH; i++){
      lattice[i] = rand() % 2;
      printf("%d,\n", lattice[i]);
    }
}
