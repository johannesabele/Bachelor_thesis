#include "heart.h"

void initiatePlanarwave(double u[size][size], double v[size][size]){
    /* sets initial conditions so that a planar wave emerges */

    //set initial values 
    for(int i=0; i < size; i++){
        for(int j=0 ; j < 7; j++){
            u[i][j] = 1;
            v[i][j] = 0;
        }
    }
}

void initiateSphericalwave(double u[size][size], double v[size][size]){
    /* sets initial conditions so that a spherical wave emerges */

    //set initial values
    for(int i=46; i<56; i++){
        for(int j=46; j<56; j++){
            u[i][j] = 1;
            v[i][j] = 0;
        }
    }
}

void initiateSpiralwave(double u[size][size], double v[size][size]){
    /* sets initial conditions for a spiral wave */

    //upper half: u=1
    for(int i=0; i < (int)(size/2); i++){
        for(int j=0; j < size; j++){
            u[i][j] = 1;
        }
    }

    //left lower quater: v=1
    for(int i=(int)(size/2); i < size; i++){
        for(int j=0; j < (int)(size/2 + 5); j++){
            v[i][j] = 1;
        }
    }
}

void initiateFourSpiralwaves(double u[size][size], double v[size][size]){

    /****************** upper half *************************/

    //upper horizontal strip: u=1
    for(int i=0; i < (int)(size/4); i++){
        for(int j=0; j < size; j++){
            u[i][j] = 1;
        }
    }

    //left lower quater: v=1
    for(int i=(int)(size/4); i < (int)(size/2); i++){
        for(int j=0; j < (int)(size/4); j++){
            v[i][j] = 1;
        }
    }
    
    //left lower quater: v=1
    for(int i=(int)(size/4); i < (int)(size/2); i++){
        for(int j=(int)(3*size/4); j < size; j++){
            v[i][j] = 1;
        }
    }

    /********************** lower half **********************/

    //lower horizontal strip: u=1
    for(int i=(int)(3*size/4); i < size; i++){
        for(int j=0; j < size; j++){
            u[i][j] = 1;
        }
    }

    //left lower vetical strip: v=1
    for(int i=(int)(size/2); i < (int)(3*size/4); i++){
        for(int j=0; j < (int)(size/4); j++){
            v[i][j] = 1;
        }
    }


    //right lower vertical strip: v=1
    for(int i= (int)(size/2); i < (int)(3*size/4); i++){
        for(int j= (int)(3*size/4); j < size; j++){
            v[i][j] = 1;
        }
    }

    return;
}

void initiateHalfShereical(double u[size][size], double v[size][size]){

    //set initial values 
    for(int i=48; i < 53; i++){
        for(int j=0 ; j < 2; j++){
            u[i][j] = 1;
            v[i][j] = 0;
        }
    }
}

void initiateLowerStrips(double u[size][size], int start){

    //lower strip
    for(int i=size-7; i < size-1; i++){
        for(int j=start-10 ; j <  start; j++){
            u[i][j] = 1;
            //v[i][j] = 0;
        }
    }
}

void initiateScrollwave(double u[size][size], double v[size][size]){
    /* sets initial conditions for a spiral wave */

    //upper half: u=1
    for(int i=0; i < (int)(size/2); i++){
        for(int j=0; j < size; j++){
            u[i][j] = 1;
        }
    }

    //left lower quater: v=1
    for(int i=(int)(size/2); i < size; i++){
        for(int j=0; j < (int)(size/2 + 30); j++){
            v[i][j] = 1;
        }
    }
}
