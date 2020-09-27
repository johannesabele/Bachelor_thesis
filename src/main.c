#include "heart.h"

/********************************* Declarations *************************************************/

void calculateDiffusion(double u[size][size], double diff[size][size]){
    /* writes diffusion into dedicated diffusion field with nine point stencil format*/

    for(int i=1; i < size-1; i++){
        for(int j=1; j < size-1; j++){
            
            //nine-point stencil
            diff[i][j] = (4*u[i+1][j] + 4*u[i-1][j] +  4*u[i][j+1] + 4*u[i][j-1] + u[i+1][j+1] + 
                         u[i+1][j-1] + u[i-1][j+1] + u[i-1][j-1] - 20*u[i][j])/6/spacing/spacing;
        }
    }
}

void enforceNoflux(double u[size][size]){
    /* enforces no-flux boundary conditions, only used in previous version */

    for(int i=0; i<size; i++){
        u[i][0] = u[i][1];
        u[i][size-1] = u[i][size-2];
    }

    for(int j=1; j<size-1; j++){
        u[0][j] = u[1][j];
        u[size-1][j] = u[size-2][j];
    }
}

void enforceDirichlet(double u[size][size]){
    /* enforces Dirichlet boundary conditions */

    for(int i=0; i<size; i++){
        u[i][0] = 0;
        u[i][size-1] = 0;
    }

    for(int j=0; j<size; j++){
        u[0][j] = 0;
        u[size-1][j] = 0;
    }
}

void alievpanfilovStep(double dt, double u[size][size], double v[size][size], 
                       double diff[size][size], int i, int j){
    /* applies explicit euler step to aliev-panfliov DGLs */

    //compute epsilon
    double epsilon = epsilon_0 + mu_1*v[i][j]/(mu_2 + u[i][j]);

    //update values at grid point i,j
    double u_current = u[i][j] + dt*(D*diff[i][j] - k*u[i][j]*(u[i][j] - a)*(u[i][j] - 1) - u[i][j]*v[i][j]);
    double v_current = v[i][j] + dt*epsilon*(-v[i][j] - k*u[i][j]*(u[i][j] - b - 1));
    u[i][j] = u_current;
    v[i][j] = v_current;
}

void activeStressStep(double dt, double u[size][size], double T_a[size][size], int i, int j){
    /* applies explicit euler step to DGL for mechnical stress */
    
    //physiologically right epsilon
    double epsilon_T = 1 - 0.9*exp(-exp(-30*(fabs(u[i][j]) - 0.1)));

    //physiologically wrong epsilon
    //double epsilon_T = 0.5 + 4.5*exp(-exp(-0.2*(u[i][j] - 0.5)));

    //update T_a field
    T_a[i][j] = T_a[i][j] + dt*epsilon_T*(k_T*fabs(u[i][j]) - T_a[i][j]);
}

void emSimulation(
    double t,
    const double N_max, 

    double u[size][size], 
    double v[size][size], 
    double diff[size][size], 
    double T[size][size],
    double dA[size][size],

    double x[size_mech][size_mech][dimension],
    double x_old[size_mech][size_mech][dimension],
    double f[size_mech][size_mech][dimension],
     
    double n[size_mech][size_mech],
    double q[size_mech][size_mech][num_q][dimension],
    double x_cm[size_mech][size_mech][dimension],
    double rest_lengths[size_mech][size_mech][2],

    double T_dy[size][size],
    double dA_dy[size][size],
    double x_dy[size_mech][size_mech][dimension],
    double x_old_dy[size_mech][size_mech][dimension],
    double f_dy[size_mech][size_mech][dimension],
    double q_dy[size_mech][size_mech][num_q][dimension],
    double x_cm_dy[size_mech][size_mech][dimension],

    FILE* electr, FILE* mech, FILE* mech_dy, FILE* par, 
    
    int part_i, int part_j,

    bool set_mech){
    /* simulate electromechanics up to N_max*/
    
    //help array
    double hlpx[size_mech][size_mech][dimension];

    //simulation 
    int cnt = 0;
    int pic_cnt = 0;
    int start = 0;
    int N_one_unit = 50;

    while(cnt < N_max){
        
        //induce chaos 
        if(chaos == 0){
            if(cnt >= 14*N_one_unit && cnt <= 30*N_one_unit && (cnt%(N_one_unit)) == 0) {
                initiateLowerStrips(u, 28 + start); 
                initiateLowerStrips(u, 28);    
                start += 4;
            }
        }

        //update electrical variables
        for(int it=0; it < it_e; it++){
            calculateDiffusion(u, diff);
            
            //update u,v,T by one time step
            for(int i=0; i < size; i++){
                for(int j=0; j < size; j++){
                    
                    //active stress first because current u is needed, after alievpanfilovStep u at next time step
                    activeStressStep(delta_t_e, u, T, i, j);
                    alievpanfilovStep(delta_t_e, u, v, diff, i, j);
                }
            }

            //enforce boundary conditions 
            enforceNoflux(u);
        }
        
        //update mechanical arrays
        if(set_mech == 0 && cnt > N_start_mech){
            
            //iterative Loop
            for(int it=0; it < it_m; it++){ 
                
                //save old postitions for next time step
                memcpy(hlpx, x, sizeof hlpx);
                
                //mechanical routine
                updateQ(q,x,n);
                updateCentersOfMass(x_cm, x);
                resetForces(f);
                updateAxialForces(f,q,x_cm,rest_lengths,n,T);
                updateStructuralForces(f, x, rest_lengths);
                updatePositions(delta_t_m, x, x_old, f);
                updateDeformation(dA, x);

                memcpy(x_old, hlpx, sizeof hlpx);
            }
        }

        //update arrays for dynamical prediction
        if(coupling == 0 && cnt > N_coupling){
            
            //update dynamical active stress with u=dA 
            for(int i=0; i < size; i++){
                for(int j=0; j < size; j++){
                    activeStressStep(delta_t_e, dA, T_dy, i, j);
                }
            }
            
            //iterative Loop
            for(int it=0; it < it_m; it++){ 
                memcpy(hlpx, x_dy, sizeof hlpx);
                
                updateQ(q_dy,x_dy,n);
                updateCentersOfMass(x_cm_dy, x_dy);
                resetForces(f_dy);
                updateAxialForces(f_dy,q_dy,x_cm_dy,rest_lengths,n,T_dy);
                updateStructuralForces(f_dy, x_dy, rest_lengths);
                updatePositions(delta_t_m, x_dy, x_old_dy, f_dy);
                updateDeformation(dA_dy, x_dy);

                memcpy(x_old_dy, hlpx, sizeof hlpx);
            }
        }

        //print to .csv files
        if(cnt >= N_output && (cnt%sample_rate) == 0){
            printArray(u, electr);
            printArray(dA, mech);
            //printArray(dA_dy, mech_dy);
            //fprintf(par, "%d, %g, %g, %g, %g\n", cnt, T[part_i][part_j], u[part_i][part_j], v[part_i][part_j], dA[part_i][part_j]);
            pic_cnt++;
            if((pic_cnt%10) == 0) printf("%d/%d\n", pic_cnt, (int)(N_max - N_output)/sample_rate);
        }
    
    cnt++;
    }
    
    return;
}


/******************** MAIN **********************/

int main(){


    /******************************************  Initialization  ******************************************/

    //particles index for printing in particle.csv
    int part_i = size/2;
    int part_j = size/4;

    //load parameters from config.ini, use "spiral" for spiral, "scroll" for scroll, "chaos" for chaos 
    char mode[] = "chaos";
    loadParams(mode);
    size_mech = size + 2*pad + 1;

    //initialize electrical arrays
    double u[size][size];
    double v[size][size];
    double T[size][size];
    double diff[size][size];
    double dA[size][size];

    //initialize mechanical arrays
    double x[size_mech][size_mech][dimension];
    double x_old[size_mech][size_mech][dimension];
    double f[size_mech][size_mech][dimension];
    double n[size_mech][size_mech];
    double q[size_mech][size_mech][num_q][dimension];
    double x_cm[size_mech][size_mech][dimension];
    double rest_lengths[size_mech][size_mech][2];

    //initialize arrays for dynamical system (_dy arrays)
    double T_dy[size][size];
    double dA_dy[size][size];
    double x_dy[size_mech][size_mech][dimension];
    double x_old_dy[size_mech][size_mech][dimension];
    double f_dy[size_mech][size_mech][dimension];
    double q_dy[size_mech][size_mech][num_q][dimension];
    double x_cm_dy[size_mech][size_mech][dimension];


    //create data file for electrical data
    FILE* electr = fopen("electrical.csv", "w");
    if(electr == NULL){
        
        printf("Error: Cannot write in file\n");
        return 0;
    }

    //create data file for mechnical data
    FILE* mech = fopen("mechanical.csv", "w");
    if(mech == NULL){
        
        printf("Error: Cannot write in file\n");
        return 0;
    }

    //create data file for dynamical data
    FILE* mech_dy = fopen("mechanical_dy.csv", "w");
    if(mech == NULL){
        
        printf("Error: Cannot write in file\n");
        return 0;
    }

    //create data file for one particle
    FILE* par = fopen("particle.csv", "w");
    if(par == NULL){
        
        printf("Error: Cannot write in file\n");
        return 0;
    }
    
    //set all to zero
    resetElectrArrays(u,v,T,diff);
    initializeMechArrays(x, x_old, f, n);

    //set dynamical arrays to zero
    resetElectrArrays(u,v,T_dy,diff);
    initializeMechArrays(x_dy, x_old_dy, f_dy, n);
    
    //set initial values
    if(strcmp(mode, "chaos") == 0)  initiatePlanarwave(u,v);
    if(strcmp(mode, "spiral") == 0)  initiateSpiralwave(u,v);
    if(strcmp(mode, "scroll") == 0)  initiateScrollwave(u,v);
    
    //preperation for mechanical grid 
    updateCentersOfMass(x_cm, x);
    updateQ(q,x,n);
    calculateRestlengths(q, x_cm, rest_lengths);
    
    //preperation for dynamical grid
    updateCentersOfMass(x_cm_dy, x_dy);
    updateQ(q_dy,x_dy,n);


    /*****************************************  Simulation  ***********************************************/

    //set time parameters
    double t = 0;
    printf("Number of pictures: %d\n", (int)((N_max - N_output)/sample_rate));
    
    //execute simulation
    emSimulation(t, N_max,
                 u, v, diff, T, dA, 
                 x, x_old, f, n, q, x_cm, rest_lengths, 
                 T_dy, dA_dy, x_dy, x_old_dy, f_dy, q_dy, x_cm_dy, 
                 electr, mech, mech_dy, par,
                 part_i, part_j,
                 set_mech);

    
    //close .csv files
    fclose(electr);
    fclose(mech);
    fclose(mech_dy);
    fclose(par);
    
    return 0;
}