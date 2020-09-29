#ifndef HEART_H
#define HEART_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

/****************************** Constants **********************************/
 
double a;                           //electrical model parameters
double b;
double mu_1;
double mu_2;
double k;
double epsilon_0;
double D;
double delta_t_e;
double it_e;

double k_T;                         //mechanical parameters
double k_ij;    
double k_ij_pad;                   
double k_j;
double k_a;
double k_a_pad;
double c_a;
double m;
double c_damp;
double delta_t_m;  
double it_m;   
 
int dimension;                      //dimension of vectors
int size;                           //size of simulation grid arrays 
int pad;                            //num padding cells
int size_mech;                      //size of mechanical grid
int num_q;                          //number of intersection points per cell
 
double n_0;                         //initial n for fiber orientation in j direction
double l_0;                         //rest length of structural springs

double spacing;                     //spacing between cells
double A_undeformed;

int N_max;                          //time parameters 
int N_start_mech;
int N_output;
int N_coupling;
int sample_rate;         

int set_mech;                       //switches for behavior
int coupling;
int chaos;



/****************************  Declaration of functions ******************************/


/*----------- initialconditions -----------*/

void initiatePlanarwave(double u[size][size], double v[size][size]);
void initiateSphericalwave(double u[size][size], double v[size][size]);
void initiateSpiralwave(double u[size][size], double v[size][size]);
void initiateFourSpiralwaves(double u[size][size], double v[size][size]);
void initiateHalfShereical(double u[size][size], double v[size][size]);
void initiateLowerStrips(double u[size][size], int start);
void initiateScrollwave(double u[size][size], double v[size][size]);

/*----------- array handling ----------------------*/

void resetElectrArrays(double u[size][size], double v[size][size], double T[size][size], 
double diff[size][size]);
void printArray(double arr[size][size], FILE* f);
void printMechArray(double arr[size_mech][size_mech][dimension], FILE *f);
void printParticle(double t, double u[size][size], double v[size][size], double T[size][size], double dA[size][size], int i, int j, FILE *f);

/*----------- vector calculus functions -----------*/

double *addVec(double vec1[dimension], double vec2[dimension]);
double *substractVec(double vec1[dimension], double vec2[dimension]);
double *multiplyVec(double c, double vec[dimension]);
double absValue(double* vec);
double *normalizeVec(double *vec);

/*------------ mechanical simulation --------------*/

//forces
double f_a(double k, double *q, double *x_cm, double l_a0, double T_a);
double f_j(double *q, double *x_cm, double l_j0, double T_a);
double f_ij(double k_ij, double *x_i, double *x_j, double l_ij);
void updateAxialForces(double f[size_mech][size_mech][dimension], double q[size_mech][size_mech][num_q][dimension], 
                       double x_cm[size_mech][size_mech][dimension], 
                       double rest_lengths[size_mech][size_mech][2], double n[size_mech][size_mech], double T[size][size]);
void updateStructuralForces(double f[size_mech][size_mech][dimension], double x[size_mech][size_mech][dimension], 
                            double rest_lengths[size_mech][size_mech][2]);

//updating arrays
void initializeMechArrays(double x[size_mech][size_mech][dimension], double x_old[size_mech][size_mech][dimension],
                          double f[size_mech][size_mech][dimension], double n[size_mech][size_mech]);
void resetForces(double f[size_mech][size_mech][dimension]);
void calculateRestlengths(double q[size_mech][size_mech][num_q][dimension], double x_cm[size_mech][size_mech][dimension], 
                          double rest_lengths[size_mech][size_mech][2]);
void updateQ(double q[size_mech][size_mech][num_q][dimension], double x[size_mech][size_mech][dimension], double n[size_mech][size_mech]);
void updateCentersOfMass(double x_cm[size_mech][size_mech][dimension], double x[size_mech][size_mech][dimension]);
void updateDisplacement(double displacement[size_mech][size_mech][dimension], 
                        double x_cm[size_mech][size_mech][dimension], double x_cm_old[size_mech][size_mech][dimension]);
void updateDeformation(double dA[size][size], double x[size][size][dimension]);

//numerical algorithm
void updatePositions(double delta_t_m, double x[size_mech][size_mech][dimension], double x_old[size_mech][size_mech][dimension], 
                    double f[size_mech][size_mech][dimension]);
 

/*--------------------- helper -----------------------*/

//load parameters 
void loadParams(char behavior[]);

#endif