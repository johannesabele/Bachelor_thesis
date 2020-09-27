#include "heart.h"
#include "ini.c"

//loads parameters from congfig.ini
void loadParams(char behaviour[]){

    printf("Loading parameters...\n");
    if(strcmp("spiral",behaviour) == 0){
        printf("Using spiral parameters...\n");
        ini_t *conf = ini_load("config.ini");
        if(conf){
            //laod eletrical parameters
            ini_sget(conf, "spiral", "a", "%lf", &a);
            ini_sget(conf, "spiral", "b", "%lf", &b);
            ini_sget(conf, "spiral", "mu_1", "%lf", &mu_1);
            ini_sget(conf, "spiral", "mu_2", "%lf", &mu_2);
            ini_sget(conf, "spiral", "k", "%lf", &k);
            ini_sget(conf, "spiral", "epsilon_0", "%lf", &epsilon_0);
            ini_sget(conf, "spiral", "D", "%lf", &D);
            ini_sget(conf, "spiral", "delta_t_e", "%lf", &delta_t_e);
            ini_sget(conf, "spiral", "it_e", "%lf", &it_e);
            //load mechanical parameters
            ini_sget(conf, "spiral", "k_T", "%lf", &k_T);
            ini_sget(conf, "spiral", "k_ij", "%lf", &k_ij);
            ini_sget(conf, "spiral", "k_ij_pad", "%lf", &k_ij_pad);
            ini_sget(conf, "spiral", "k_j", "%lf", &k_j);
            ini_sget(conf, "spiral", "k_a", "%lf", &k_a);
            ini_sget(conf, "spiral", "k_a_pad", "%lf", &k_a_pad);
            ini_sget(conf, "spiral", "c_a", "%lf", &c_a);
            ini_sget(conf, "spiral", "c_damp", "%lf", &c_damp);
            ini_sget(conf, "spiral", "delta_t_m", "%lf", &delta_t_m);
            ini_sget(conf, "spiral", "it_m", "%lf", &it_m);
            //load sizes of arrays
            ini_sget(conf, "spiral", "dimension", "%d", &dimension);
            ini_sget(conf, "spiral", "size", "%d", &size);
            ini_sget(conf, "spiral", "pad", "%d", &pad);
            ini_sget(conf, "spiral", "size_mech", "%d", &size_mech);
            ini_sget(conf, "spiral", "num_q", "%d", &num_q);
            //load fibre orientation
            ini_sget(conf, "spiral", "n_0", "%lf", &n_0);
            ini_sget(conf, "spiral", "l_0", "%lf", &l_0);
            ini_sget(conf, "spiral", "spacing", "%lf", &spacing);
            ini_sget(conf, "spiral", "A_undeformed", "%lf", &A_undeformed);
            //load time parameters
            ini_sget(conf, "spiral", "N_max", "%d", &N_max);
            ini_sget(conf, "spiral", "N_start_mech", "%d", &N_start_mech);
            ini_sget(conf, "spiral", "N_coupling", "%d", &N_coupling);
            ini_sget(conf, "spiral", "N_output", "%d", &N_output);
            ini_sget(conf, "spiral", "sample_rate", "%d", &sample_rate);
            //load parameters for behaviour
            ini_sget(conf, "spiral", "set_mech", "%d", &set_mech);
            ini_sget(conf, "spiral", "coupling", "%d", &coupling);
            ini_sget(conf, "spiral", "chaos", "%d", &chaos);
        }
        else{
            printf("Error: could not read config.ini\n");
        }
    }

    else if(strcmp("scroll",behaviour) == 0){
        printf("Using scroll parameters...\n");
        ini_t *conf = ini_load("config.ini");
        if(conf){
            //laod eletrical parameters
            ini_sget(conf, "scroll", "a", "%lf", &a);
            ini_sget(conf, "scroll", "b", "%lf", &b);
            ini_sget(conf, "scroll", "mu_1", "%lf", &mu_1);
            ini_sget(conf, "scroll", "mu_2", "%lf", &mu_2);
            ini_sget(conf, "scroll", "k", "%lf", &k);
            ini_sget(conf, "scroll", "epsilon_0", "%lf", &epsilon_0);
            ini_sget(conf, "scroll", "D", "%lf", &D);
            ini_sget(conf, "scroll", "delta_t_e", "%lf", &delta_t_e);
            ini_sget(conf, "scroll", "it_e", "%lf", &it_e);
            //load mechanical parameters
            ini_sget(conf, "scroll", "k_T", "%lf", &k_T);
            ini_sget(conf, "scroll", "k_ij", "%lf", &k_ij);
            ini_sget(conf, "scroll", "k_ij_pad", "%lf", &k_ij_pad);
            ini_sget(conf, "scroll", "k_j", "%lf", &k_j);
            ini_sget(conf, "scroll", "k_a", "%lf", &k_a);
            ini_sget(conf, "scroll", "k_a_pad", "%lf", &k_a_pad);
            ini_sget(conf, "scroll", "c_a", "%lf", &c_a);
            ini_sget(conf, "scroll", "c_damp", "%lf", &c_damp);
            ini_sget(conf, "scroll", "delta_t_m", "%lf", &delta_t_m);
            ini_sget(conf, "scroll", "it_m", "%lf", &it_m);
            //load sizes of arrays
            ini_sget(conf, "scroll", "dimension", "%d", &dimension);
            ini_sget(conf, "scroll", "size", "%d", &size);
            ini_sget(conf, "scroll", "pad", "%d", &pad);
            ini_sget(conf, "scroll", "size_mech", "%d", &size_mech);
            ini_sget(conf, "scroll", "num_q", "%d", &num_q);
            //load fibre orientation
            ini_sget(conf, "scroll", "n_0", "%lf", &n_0);
            ini_sget(conf, "scroll", "l_0", "%lf", &l_0);
            ini_sget(conf, "scroll", "spacing", "%lf", &spacing);
            ini_sget(conf, "scroll", "A_undeformed", "%lf", &A_undeformed);
            //load time parameters
            ini_sget(conf, "scroll", "N_max", "%d", &N_max);
            ini_sget(conf, "scroll", "N_start_mech", "%d", &N_start_mech);
            ini_sget(conf, "scroll", "N_coupling", "%d", &N_coupling);
            ini_sget(conf, "scroll", "N_output", "%d", &N_output);
            ini_sget(conf, "scroll", "sample_rate", "%d", &sample_rate);
            //load parameters for behaviour
            ini_sget(conf, "scroll", "set_mech", "%d", &set_mech);
            ini_sget(conf, "scroll", "coupling", "%d", &coupling);
            ini_sget(conf, "scroll", "chaos", "%d", &chaos);
        }
        else{
            printf("Error: could not read config.ini\n");
        }
    }

    else if(strcmp("chaos",behaviour) == 0){
        ini_t *conf = ini_load("config.ini");
        printf("Using chaos parameters...\n");
        if(conf){
            //laod eletrical parameters
            ini_sget(conf, "chaos", "a", "%lf", &a);
            ini_sget(conf, "chaos", "b", "%lf", &b);
            ini_sget(conf, "chaos", "mu_1", "%lf", &mu_1);
            ini_sget(conf, "chaos", "mu_2", "%lf", &mu_2);
            ini_sget(conf, "chaos", "k", "%lf", &k);
            ini_sget(conf, "chaos", "epsilon_0", "%lf", &epsilon_0);
            ini_sget(conf, "chaos", "D", "%lf", &D);
            ini_sget(conf, "chaos", "delta_t_e", "%lf", &delta_t_e);
            ini_sget(conf, "chaos", "it_e", "%lf", &it_e);
            //load mechanical parameters
            ini_sget(conf, "chaos", "k_T", "%lf", &k_T);
            ini_sget(conf, "chaos", "k_ij", "%lf", &k_ij);
            ini_sget(conf, "chaos", "k_ij_pad", "%lf", &k_ij_pad);
            ini_sget(conf, "chaos", "k_j", "%lf", &k_j);
            ini_sget(conf, "chaos", "k_a", "%lf", &k_a);
            ini_sget(conf, "chaos", "k_a_pad", "%lf", &k_a_pad);
            ini_sget(conf, "chaos", "c_a", "%lf", &c_a);
            ini_sget(conf, "chaos", "c_damp", "%lf", &c_damp);
            ini_sget(conf, "chaos", "delta_t_m", "%lf", &delta_t_m);
            ini_sget(conf, "chaos", "it_m", "%lf", &it_m);
            //load sizes of arrays
            ini_sget(conf, "chaos", "dimension", "%d", &dimension);
            ini_sget(conf, "chaos", "size", "%d", &size);
            ini_sget(conf, "chaos", "pad", "%d", &pad);
            ini_sget(conf, "chaos", "size_mech", "%d", &size_mech);
            ini_sget(conf, "chaos", "num_q", "%d", &num_q);
            //load fibre orientation
            ini_sget(conf, "chaos", "n_0", "%lf", &n_0);
            ini_sget(conf, "chaos", "l_0", "%lf", &l_0);
            ini_sget(conf, "chaos", "spacing", "%lf", &spacing);
            ini_sget(conf, "chaos", "A_undeformed", "%lf", &A_undeformed);
            //load time parameters
            ini_sget(conf, "chaos", "N_max", "%d", &N_max);
            ini_sget(conf, "chaos", "N_start_mech", "%d", &N_start_mech);
            ini_sget(conf, "chaos", "N_coupling", "%d", &N_coupling);
            ini_sget(conf, "chaos", "N_output", "%d", &N_output);
            ini_sget(conf, "chaos", "sample_rate", "%d", &sample_rate);
            //load parameters for behaviour
            ini_sget(conf, "chaos", "set_mech", "%d", &set_mech);
            ini_sget(conf, "chaos", "coupling", "%d", &coupling);
            ini_sget(conf, "chaos", "chaos", "%d", &chaos);
        }
        else{
            printf("Error: could not read config.ini\n");
        }
    }

    else {
        printf("No parameters found! Change behaviour!\n");
        return;
    }

    printf("All parameters loaded!\n");
}