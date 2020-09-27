#include "heart.h"

double* addVec(double vec1[dimension], double vec2[dimension]){
    /* adds two 2D arrays and returns result */

    static double result[2];

    result[0] = vec1[0] + vec2[0];
    result[1] = vec1[1] + vec2[1];

    return result;
}

double* substractVec(double vec1[dimension], double vec2[dimension]){
    /* substracts vec2 from vec1 */

    static double result[2];

    result[0] = vec1[0] - vec2[0];
    result[1] = vec1[1] - vec2[1];

    return result;
}

double* multiplyVec(double c, double vec[dimension]){
    /* multiplies each entry with constant */

    static double result[2];
    
    result[0] = c*vec[0];
    result[1] = c*vec[1];

    return result;
}

double absValue(double* vec){
    /* norm in 2D */

    return sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
}

double *normalizeVec(double *vec){
    /* returns normalized vector as 2d array */

    static double result[2];

    double c = absValue(vec);

    if( c < 0.00000000000000001){
        printf("ERROR: Division through 0!\n");
        return result;
    }

    result[0] = vec[0]/c;
    result[1] = vec[1]/c;

    return result; 
}
