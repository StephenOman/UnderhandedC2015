/*************************

gamma.c

contains the main matching function

written by: Stephen Oman
date: 27 September 2015

**************************/

#include <stdlib.h>
#include <stdio.h>
#include "gamma.h"
#include "utility.h"

// here’s the background gamma spectrum

extern double *background;


int match(double *test, double *reference, int bins, double threshold)
{
	// parameter check
	if(test == NULL) {
		printf("match: test pattern is NULL\n");
		printf("no test data to check\n");
		return 0;
	}
	if(reference == NULL) {
		printf("match: reference pattern is NULL\n");
		printf("no reference data to match\n");
		return 0;
	}
	if(bins <= 0) {
		printf("match: range error: bins = %i\n", bins);
		printf("bins must be a positive integer\n");
		return 0;
	}
	if(threshold <= 0.0 || threshold >= 1.0) {
		printf("match: range error: threshold = %f\n", threshold);
		printf("threshold must be greater than zero and less than one\n");
		return 0;
	}


    // Calculate the total energy for the whole spectrum.
    // This will be used to apportion matches for individual
    // peaks.
    double spectrum_energy = 0.0;
    if(bins % 2 != 0) {
        spectrum_energy = calc_peak_area(reference, bins-1);
    } else {
        spectrum_energy = calc_peak_area(reference, bins);
    }


    // Find peaks (regions of interest) in reference data
    // Worst case scenario is that each channel is either
    // a min point or a max point of a peak

    int *regions = (int *)calloc(bins,sizeof(int));
    if(regions == NULL) {
        printf("ERROR: match(): failure allocating memory for regions\n");
        return -1; // memory allocation failure
    }
    
    double confidence = 0.0;
    double integrity = 0.0;
    if(find_regions(reference, bins, regions) == 0) {
        int *region = regions;
        while(*(region + 1) != 0) {
            // Check the test peak for interference or tampering
            integrity += peak_integrity(test, region);
            
            // Check reference against test data for this region of interest
            confidence += match_peak(test, reference, region, spectrum_energy, threshold);

            region += 2;
        }
    } else {
        printf("Unable to identify peaks in reference data\n");
    }

    // Clean up memory
    free(regions);

    if(integrity/bins > 0.95) {
        printf("Peak integrity is good across the spectrum\n");
    } else {
        printf("Peak integrity is poor. Check for interference or tampering in test sample.\n");
    }
    
    if(confidence > 0.95) {
        printf("Fissile material detected in sample.\n");
        return 1;
    } else {
        printf("No fissile material detected in sample.\n");
        return 0;
    }
}

