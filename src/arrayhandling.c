#include "heart.h"

void resetElectrArrays(double u[size][size], double v[size][size], double T[size][size], 
                        double diff[size][size]){
    /* initialize u,v and T arrays */

    //set all entries in u, v and T to zero
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            
            u[i][j] = 0;
            v[i][j] = 0;
            diff[i][j] = 0;
            T[i][j] = 0;
        }
    }
}


void printArray(double arr[size][size], FILE* f){
    /* prints complete array in file (with boundary conditions) */

    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            if(j != size-1){
                fprintf(f, "%g,", arr[i][j]);
            }
            else{
                fprintf(f, "%g", arr[i][j]);
            }
        }
    fprintf(f, "\n");
    }

    fprintf(f, "\n");
}


void printMechArray(double arr[size_mech][size_mech][dimension], FILE *f){
    /* prints positions of mechanical array to file */

    for(int i=pad; i < (size_mech - pad); i++){
        for(int k=0; k < dimension; k++){
            for(int j=pad; j < (size_mech - pad); j++){            
                fprintf(f, "%g,", arr[i][j][k]);
            }
        fprintf(f, "\n");
        }
    fprintf(f, "\n\n");
    }

    fprintf(f, "\n\n\n\n");
}

void printParticle(double t, double u[size][size], double v[size][size], double T[size][size], double dA[size][size], int i, int j, FILE *f){
    /* prints variables in file in format t u v T dA */

    fprintf(f, "%g %g %g %g %g\n", t, u[i][j], v[i][j], T[i][j], dA[i][j]);
}