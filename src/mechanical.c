#include "heart.h"


/* ----------------------- Forces ------------------------------*/

double f_a(double k, double *q, double *x_cm, double l_a0, double T_a){
    /* calculates absolute value of force from active spring */

    return -k*(absValue(substractVec(q,x_cm)) - l_a0/(1+ c_a*T_a));
}

double f_j(double *q, double *x_cm, double l_j0, double T_a){
    /* calculates absolute value of force from non-active spring */

    return -k_j*(absValue(substractVec(q, x_cm)) - l_j0);
}

double f_ij(double k, double *x_i, double *x_j, double l_ij){
    /* calculates absolute value of structural forces */

    return -k*(absValue(substractVec(x_i,x_j)) - l_ij);
}

void updateAxialForces(
    double f[size_mech][size_mech][dimension], 
    double q[size_mech][size_mech][num_q][dimension], 
    double x_cm[size_mech][size_mech][dimension], 
    double rest_lengths[size_mech][size_mech][2], 
    double n[size_mech][size_mech], 
    double T[size][size]){
    /* updates forces from axial springs */

    double local_force;
    double local_T;
    double local_k_a;

    for(int i=0; i < size_mech-1; i++){
        for(int j=0; j < size_mech-1; j++){
            
            local_T = 0;
            local_k_a = k_a_pad;

            //enforce same T for two rows of cells above and below simulation grid
            if(i == (pad-1) && (j >= pad) && (j < size+pad)){
                local_T = T[i-pad+1][j-pad];
                local_k_a = k_a;
            }
            if(i == (pad-2) && (j >= pad) && (j < size+pad)){
                local_T = T[i-pad+2][j-pad];
                local_k_a = k_a;
            }
            if(i == (size + pad +1) && (j >= pad) && (j < size+pad)){
                local_T = T[i-pad-2][j-pad];
                local_k_a = k_a;
            }
            if(i == (size+pad) && (j >= pad) && (j < size+pad)){
                local_T = T[i-pad-1][j-pad];
                local_k_a = k_a;
            }

            //enforce same T for two rows of cells right and left of simulation grid
            if(j == (pad-1) && (i >= pad) && (i < size+pad)){
                local_T = T[i-pad][j-pad+1];
                local_k_a = k_a;
            }
            if(j == (pad-2) && (i >= pad) && (i < size+pad)){
                local_T = T[i-pad][j-pad+2];
                local_k_a = k_a;
            }
            if(j == (size + pad +1) && (i >= pad) && (i < size+pad)){
                local_T = T[i-pad][j-pad-2];
                local_k_a = k_a;
            }
            if(j == (size+pad) && (i >= pad) && (i < size+pad)){
                local_T = T[i-pad][j-pad-1];
                local_k_a = k_a;
            }
            
            //enfore different k_a, T_a for simulation grid
            if((i >= pad) && (i < size+pad) && (j >= pad) && (j < size+pad)){
                local_T = T[i-pad][j-pad];                     
                local_k_a = k_a;
            }
            
            for(int k=0; k < dimension; k++){
                
                //distribute force from q_0 spring 
                local_force = f_a(local_k_a, q[i][j][0], x_cm[i][j], rest_lengths[i][j][0], local_T)*normalizeVec(substractVec(q[i][j][0], x_cm[i][j]))[k];
                f[i][j][k] += (1-n[i][j])*local_force;
                f[i+1][j][k] += n[i][j]*local_force;

                //distribute force from q_1
                local_force = f_j(q[i][j][1], x_cm[i][j], rest_lengths[i][j][1], local_T)*normalizeVec(substractVec(q[i][j][1], x_cm[i][j]))[k];
                f[i][j][k] += n[i][j]*local_force;
                f[i][j+1][k] += (1-n[i][j])*local_force;

                //distribute force from q_2
                local_force = f_a(local_k_a, q[i][j][2], x_cm[i][j], rest_lengths[i][j][0], local_T)*normalizeVec(substractVec(q[i][j][2], x_cm[i][j]))[k];
                f[i][j+1][k] += n[i][j]*local_force;
                f[i+1][j+1][k] += (1-n[i][j])*local_force;

                //distribute force from q_3
                local_force = f_j(q[i][j][3], x_cm[i][j], rest_lengths[i][j][1], local_T)*normalizeVec(substractVec(q[i][j][3], x_cm[i][j]))[k];
                f[i+1][j+1][k] += n[i][j]*local_force;
                f[i+1][j][k] += (1-n[i][j])*local_force;
            }
        }
    }    
}

void updateStructuralForces(
    double f[size_mech][size_mech][dimension], 
    double x[size_mech][size_mech][dimension], 
    double rest_lengths[size_mech][size_mech][2]){
    /* updates forces from structural springs */

    double local_force;
    double local_k;
        
    for(int i=0; i < size_mech-1; i++){
        for(int j=0; j < size_mech-1; j++){
            for(int k=0; k < dimension; k++){
                
                local_k = k_ij_pad;

                if((i >= pad) && (i < size+pad) && (j >= pad) && (j < size+pad)){
                    local_k = k_ij;
                } 

                //update force between x_0 and x_1
                local_force = f_ij(local_k, x[i][j], x[i][j+1], l_0)*normalizeVec(substractVec(x[i][j+1], x[i][j]))[k];
                f[i][j+1][k] += local_force;
                f[i][j][k] -= local_force;
            

                //update force between x_3 and x_0
                local_force = f_ij(local_k, x[i][j], x[i+1][j], l_0)*normalizeVec(substractVec(x[i+1][j], x[i][j]))[k];
                f[i+1][j][k] += local_force;
                f[i][j][k] -= local_force;
                
                
                if(j == (size_mech-2)){
                    //update force between x_1 and x_2
                    local_force = f_ij(local_k, x[i][j+1], x[i+1][j+1], l_0)*normalizeVec(substractVec(x[i+1][j+1], x[i][j+1]))[k];
                    f[i+1][j+1][k] += local_force;
                    f[i][j+1][k] -= local_force;
                }

                if(i == (size_mech-2)){
                    //update force between x_2 and x_3
                    local_force = f_ij(local_k, x[i+1][j], x[i+1][j+1], l_0)*normalizeVec(substractVec(x[i+1][j+1], x[i+1][j]))[k];
                    f[i+1][j][k] -= local_force;
                    f[i+1][j+1][k] += local_force;
                }
            }
        }
    }
}


/* --------------------------- Simulation Functions ------------------------------ */ 

void initializeMechArrays(
    double x[size_mech][size_mech][dimension],
    double x_old[size_mech][size_mech][dimension],
    double f[size_mech][size_mech][dimension],
    double n[size_mech][size_mech]){
    /* initialize size_mech arrays */ 
    
    for(int i=0; i < size_mech; i++){
        for(int j=0; j < size_mech; j++){
                
                //set positions
                x[i][j][0] = (double)j*spacing;
                x[i][j][1] = (double)i*spacing;
                x_old[i][j][0] = (double)j*spacing;
                x_old[i][j][1] = (double)i*spacing;
                
                //set one n for all cells
                n[i][j] = n_0;
                
                //reset forces
                f[i][j][0] = 0;
                f[i][j][1] = 0;
        }
    }
}

void resetForces(double f[size_mech][size_mech][dimension]){
    /* sets all entries of force array to zero */

    for(int i=0; i < size_mech; i++){
        for(int j=0; j < size_mech; j++){
            f[i][j][0] = 0;
            f[i][j][1] = 0;
        }
    }
}

void updateQ(double q[size_mech][size_mech][num_q][dimension], double x[size_mech][size_mech][dimension], double n[size_mech][size_mech]){
    /* updates the q-array */

    for(int i=0; i < size_mech-1; i++){
        for(int j=0; j < size_mech-1; j++){
            for(int k=0; k < dimension; k++){
                q[i][j][0][k] = addVec(x[i][j], multiplyVec(n[i][j], substractVec(x[i+1][j],x[i][j])))[k];
                q[i][j][1][k] = addVec(x[i][j], multiplyVec((1-n[i][j]), substractVec(x[i][j+1],x[i][j])))[k];
                q[i][j][2][k] = addVec(x[i][j+1], multiplyVec((1-n[i][j]), substractVec(x[i+1][j+1],x[i][j+1])))[k];
                q[i][j][3][k] = addVec(x[i+1][j], multiplyVec(n[i][j], substractVec(x[i+1][j+1],x[i+1][j])))[k];
            }
        }
    }
}

void calculateRestlengths(double q[size_mech][size_mech][num_q][dimension], double x_cm[size_mech][size_mech][dimension], 
    double rest_lengths[size_mech][size_mech][2]){
    /* writes rest lengths in arrays, index 0 for active spring, 1 for non-active spring */

    for(int i=0; i < size_mech-1; i++){
        for(int j=0; j < size_mech-1; j++){

            //length of active spring 
            rest_lengths[i][j][0] = absValue(substractVec(q[i][j][0], x_cm[i][j]));

            //length of non active spring 
            rest_lengths[i][j][1] = absValue(substractVec(q[i][j][1], x_cm[i][j]));
        }
    }
}

void updateCentersOfMass(double x_cm[size_mech][size_mech][dimension], double x[size_mech][size_mech][dimension]){
    /* updates the positions of the centres of mass */

    for(int i=0; i < size_mech-1; i++){
        for(int j=0; j < size_mech-1; j++){
            for(int k=0; k < dimension; k++){
                x_cm[i][j][k] = (x[i][j][k] + x[i][j+1][k] + x[i+1][j][k] + x[i+1][j+1][k])/4;
            }
        }
    }
}

void updateDisplacement(double displacement[size_mech][size_mech][dimension], double x_cm[size_mech][size_mech][dimension], double x_cm_old[size_mech][size_mech][dimension]){
    /*updates the Displacement vector field */

    for(int i=0; i < size_mech-1; i++){
        for(int j=0; j < size_mech-1; j++){
            for (int k = 0; k < dimension; k++){
                displacement[i][j][k] = substractVec(x_cm[i][j], x_cm_old[i][j])[k];                
            }
        }
    }
}

void updatePositions(
    double dt, 
    double x[size_mech][size_mech][dimension], 
    double x_old[size_mech][size_mech][dimension],
    double f[size_mech][size_mech][dimension]){
    /* updates positions with verlet alghorithm */

    for(int i=1; i < size_mech-1; i++){
        for(int j=1; j < size_mech-1; j++){
            for(int k=0; k < dimension; k++){
                x[i][j][k] = (2*x[i][j][k] - x_old[i][j][k]*(1 - c_damp*dt/2) + f[i][j][k]*dt*dt)/(1 + c_damp*dt/2);        
            }
        }
    }  
}

void updateDeformation(double dA[size][size], double x[size_mech][size_mech][dimension]){
    /* updates the deformation array with formula for generic quadrilateral */

    double a;
    double b;
    double c;
    double d;
    double diagonal1;
    double diagonal2;
    double hlp;

    for(int i=pad; i < (size + pad); i++){
        for(int j=pad; j < (size + pad); j++){
            
            dA[i-pad][j-pad] = 0;

            a = absValue(substractVec(x[i+1][j+1], x[i+1][j]));
            b = absValue(substractVec(x[i+1][j+1], x[i][j+1]));
            c = absValue(substractVec(x[i][j+1], x[i][j]));
            d = absValue(substractVec(x[i+1][j], x[i][j]));
            diagonal1 = absValue(substractVec(x[i][j+1], x[i+1][j]));
            diagonal2 = absValue(substractVec(x[i+1][j+1], x[i][j]));
            hlp = (b*b + d*d - a*a - c*c);

            dA[i-pad][j-pad] = sqrt(4*diagonal1*diagonal1*diagonal2*diagonal2 - hlp*hlp)/4;
            dA[i-pad][j-pad] = dA[i-pad][j-pad]/A_undeformed - 1;
            
            //error if calculated area is too small
            if(dA[i-pad][j-pad] < -0.998) printf("Error: Cell with almost no area\n");
        }
    }
}

