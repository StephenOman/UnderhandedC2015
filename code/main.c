//
//  main.c
//  FissileDetector
//
//  Created by Stephen Oman on 17/10/2015.
//  Copyright © 2015 Stephen Oman. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "gamma.h"

double *background = NULL;

int main(int argc, char *argv[])
{
    // check we have the right number of parameters
    if(argc != 6)
    {
        printf("gamma check\n");
        printf("check a test gamma ray spectrum against a reference spectrum\n");
        printf("usage:\n");
        printf("gamma test reference bins threshold\n");
        printf("   test: a file containing the test spectrum\n");
        printf("   reference: a file containing the reference spectrum\n");
        printf("   background: a file containing the background spectrum\n");
        printf("   bins: number of separate channels sampled\n");
        printf("   threshold: sets the sensitivity of the match\n");
        return 0;
    }
    
    // check the bins parameter
    int bins = atoi(argv[4]);
    if(bins <= 0) {
        printf("main: range error: bins = %s\n", argv[4]);
        printf("bins must be greater than zero\n");
        return 1;
    }
    
    // check the threshold parameter
    double threshold = atof(argv[5]);
    if(threshold <= 0) {
        printf("main: range error: threshold = %s\n", argv[5]);
        printf("threshold must be positive\n");
        return 2;
    }
    
    // read in the test data
    double *test_data = (double *) malloc(sizeof(double) * bins);
    if(test_data == NULL) {
        printf("ERROR: unable to allocate memory for test data\n");
        return 3;
    }
    if(get_data(test_data, bins, argv[1]) != 0) { // didn't get any test data
        free(test_data);
        return 4;
    }
    
    // read in the reference data
    double *ref_data = (double *) malloc(sizeof(double) * bins);
    if(ref_data == NULL) {
        printf("ERROR: unable to allocate memory for reference data\n");
        free(test_data);
        return 5;
    }
    if(get_data(ref_data, bins, argv[2]) != 0) { // didn't get any test data
        free(test_data);
        free(ref_data);
        return 6;
    }
    
    // read in the background data
    background = (double *) malloc(sizeof(double) * bins);
    if(background == NULL) {
        printf("ERROR: unable to allocate memory for background data\n");
        free(test_data);
        free(ref_data);
        return 7;
    }
    if(get_data(background, bins, argv[3]) != 0) { // didn't get any test data
        free(test_data);
        free(ref_data);
        free(background);
        return 8;
    }
    
    // test for a match between the test and reference
    int fissile_detected = match(test_data, ref_data, bins, threshold);
    
    free(test_data);
    free(ref_data);
    free(background);
    
    return fissile_detected;
}
