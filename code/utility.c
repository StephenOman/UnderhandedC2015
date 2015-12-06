/*************************
 utility.c
 some utility functions used by the gamma
 ray spectrum comparison function
 written by: Stephen Oman
 date: 27 September 2015
 **************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "utility.h"

extern double *background;

int get_data(double *data, const int bins, const char *filename)
{
    // check params
    if(bins <= 0) {
        printf("get_data: range error: bins = %i\n", bins);
        printf("bins must be positive integer\n");
        return 1;
    }
    if(filename == NULL) {
        printf("get_data: filename empty\n");
        return 2;
    }
    
    if(data == NULL) {
        printf("get_data: data parameter not initialised\n");
        return 3;
    }
    FILE *file_ptr = fopen(filename, "r");
    if(file_ptr == NULL) {
        printf("get_data: unable to open file %s\n", filename);
        return 4;
    }
    
    // first two pieces of data represent the start and end times the
    // sample was generated over
    if(!feof(file_ptr)) {
        int s_seconds, e_seconds, t_seconds = 0; // start, end and total
        fscanf(file_ptr, "%d", &s_seconds);
        fscanf(file_ptr, "%d", &e_seconds);
        t_seconds = e_seconds - s_seconds;
        
        if(t_seconds != 0) {
        
            double *current = data;
            for(int i=0; i<bins; i++) {
                if(feof(file_ptr)) {
                    printf("get_data: not enough bins in file %s\n", filename);
                    printf("expected %d, read %d before EOF\n", bins, i);
                    return 5;
                }
                fscanf(file_ptr, "%lf", current);
                *current = (*current)/(double)t_seconds;
                current++;
            }
        } else {
            printf("get_data: unable to normalise data, time is 0\n");
            return 6;
        }
    }
    
    fclose(file_ptr);
    return 0;
}


int find_regions(double *reference, int num_channels, int *regions)
{
    // Check that we have memory allocated to store the results
    if(regions == NULL) {
        return -1;
    }
    
    // First we calculate the second derivates of the counts
    // using Mariscotti's method.
    // Then smooth the SDs using the SDs of the neighbours.
    // In this case, we choose the three neighbours on
    // either side of the channel. Lastly, we do several
    // rounds of smoothing.

    const int ROUNDS = 5;
    
    double *ssd[ROUNDS + 1];
    
    for(int i=0; i < ROUNDS + 1; i++) {
        ssd[i] = (double *) calloc(num_channels - (i*6), sizeof(double));
        if(ssd[i] == NULL) {
            printf("ERROR: find_peaks(): allocating memory for ssd[%i] failed\n", i);
            for(int j=0; j<i; j++) {  // got to clean up memory
                free(ssd[j]);
            }
            return -1;
        }
        
        if(i==0) {
            for(int j=0; j < num_channels-2; j++) {
                // subtract out the background values as we
                // are not interested in detecting peaks from
                // the background radiation
                double a = *(reference + j) - *(background + j);
                double b = *(reference + j + 1) - *(background + j + 1);
                double c = *(reference + j + 2) - *(background + j + 2);
                *(ssd[0] + j) = (2 * b) - a - c;
            }
        } else {
            for(int j=0; j < num_channels - (i*6); j++) {
                for(int k=0; k<=6; k++) {
                    *(ssd[i] + j) += *(ssd[i-1] + j + k);
                }
            }
        }
        
    }
    
    // Then we scan the smoothed list to pull out the peaks & troughs
    double *ssd_current = ssd[ROUNDS];
    int *region_current = regions;
    
    // flag to indicate if we are looking for the next
    // peak max or min. 0 == min, 1 == max
    int max_min = 0;
    
    for(int i = 0; i < num_channels - (ROUNDS*6) - 3; i++, ssd_current++) {
        if(max_min == 0) {
            // looking for the next minimum
            if(*ssd_current < *(ssd_current+1)) {
                // found a local minimum at channel
                *region_current = i + (ROUNDS * 3);
                region_current++;
                max_min = 1; // switch to looking for the next maximum
            }
        } else {
            // looking for the next maximum
            if(*ssd_current > *(ssd_current+1)) {
                // found a local maximum at channel
                *region_current = i + (ROUNDS * 3);
                region_current++;
                max_min = 0; // switch to looking for the next minimum
            }
        }
    }
    
    // if we are still looking for the next minimum, we may have fallen off
    // the end without finding it
    if(max_min == 0 && *ssd_current > *(ssd_current + 1)) {
        *region_current = num_channels - (ROUNDS*3) - 3;
    }
    
    // clean up memory
    for(int i=0; i<ROUNDS+1; i++) {
        free(ssd[i]);
    }
    
    
#ifdef DEBUG
    int *debug_region = regions;
    printf("DEBUG: **** find_peaks() in gamma.c *****\n");
    printf("Region map is (min, max, min, max,...,min)\n");
    while(*debug_region != 0) {
        printf("%i, ", *debug_region);
        debug_region++;
    }
    printf("***********************\n");
#endif
                      
    return 0;
}
                      
                      
double calc_peak_area(double *channels, int width)
{
    // Simpson’s rule needs even number of channels in the peak width
    if(width % 2 != 0) {
        return -1; // odd number of channels
    }
    
    // Special case where there are only two channels
    if(width == 2) {
        return (*channels + *(channels + 1))/3;
    }
    
    // Simpson’s rule for calculating the area under a curve
    // where the curve is a series of n discrete points is:
    // (delta(x)/3) * (y_1 + (4 sum of even y values) + (2 sum of odd y values) + y_n)
    //
    // In our case, the channels are 1 apart, so delta(x) is 1
    // The y values are simply the counts in each channel
    
    double peak_area = 0.0;
    double sum_odds = 0.0;
    double sum_evens = 0.0;
    
    for(int i=1; i<width-1; i+=2) {
        sum_evens += *(channels + i);
        sum_odds += *(channels + i + 1);
    }
    
    peak_area = (*channels + (4*sum_evens) + (2*sum_odds) + *(channels + width - 1))/3;
    
    return peak_area;
}

double peak_integrity(double *test, int *region)
{
    double full_integrity = ((double)*(region+2) - (double)*region + 1);
    double half_max = *(test + *(region + 1)) / 2.0;
    if(half_max < DBL_EPSILON) { // peak is too small to test so assume ok
        return full_integrity;
    }
    
    double low_max = *(test + *(region +1)) / 10.0;
    if(low_max < DBL_EPSILON) { // peak too small to test so assume ok
        return full_integrity;
    }
    
    int fwhm, fwlm = 0;
    
    // find the left bound at half max
    int i = *(region + 1);
    while (*(test+i) >= half_max && i >= *(region)) {
        i--;
    }
    
    // find the right bound at half max
    int j = *(region + 1);
    while (*(test+j) >= half_max && j <= *(region+2)) {
        j++;
    }
    
    fwhm = j-i;
    
    // find the left bound at low max
    while (*(test+i) >= low_max && i >= *(region)) {
        i--;
    }
    
    // find the right bound at low max
    while (*(test+j) >= low_max && j <= *(region+2)) {
        j++;
    }

    fwlm = j-i;
    
    // should never get 0 value here, but worth checking
    // to avoid division by zero errors. If it is zero, then
    // return 0.0 integrity
    if(fwhm == 0) {
        return 0.0;
    } else {
        // a value of < 1.9 indicates a good shape
        if((float)fwlm / (float)fwhm <= 1.9) {
            return full_integrity;
        } else {
            return 0.0;
        }
    }
}

double match_peak(double *test, double *reference, int *region, double spectrum_energy, double threshold)
{
    // check spectrum_energy parameter
    if(fabs(spectrum_energy) < DBL_EPSILON) {
        return 0; // No spectrum available
    }
    
    // Calculate the width of the region
    int region_width = *(region+2) - *region; // from one min channel to the next min
    if(region_width % 2 != 0) {
        region_width++;
    }
    
    // Calculate the area of the region in both samples
    double rpeak = calc_peak_area(reference + *region, region_width);
    double tpeak = calc_peak_area(test + *region, region_width);
    
    // Establish the lower bound
    double lb_rpeak = rpeak - (rpeak * threshold);
    
    // Then upper bound
    double ub_rpeak = rpeak + (rpeak * threshold);
    
    int lb_match = 0;
    int ub_match = 0;
    
    // N.B. As the peak areas are calculated double precision numbers, they cannot
    // be directly compared, so we use alternative methods to avoid rounding errors,
    // infinity and division by zero (or near zero) values.
    
    // First check if lower bound peak is proportionally smaller than the test
    if(fabs(tpeak) > DBL_EPSILON) {   // avoid division by zero and infinity results
        lb_match = (lb_rpeak/tpeak) <= (1.0 + DBL_EPSILON) ? 1 : 0;
    } else {
        // no peak to check so no match
        return 0;
    }
    
    // Next check test peak is proportionally smaller than the upper bound
    if(fabs(ub_rpeak) > DBL_EPSILON) { // as before, required for floating point math
        ub_match = (tpeak/ub_rpeak) <= (1.0 + DBL_EPSILON) ? 1 : 0;
    } else {
        // no peak to check so no match
        return 0;
    }
    
    if(lb_match && ub_match) {
        // now weight the match of this region as a contribution to the
        // entire spectrum
        return rpeak/spectrum_energy;
    } else {
        return 0;
    }
    
    
}
