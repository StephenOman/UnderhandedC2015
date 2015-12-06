/*************************

utility.h

utility functions for the gamma ray spectroscopy check

written by: Stephen Oman
date: 27 September 2015

**************************/

#ifndef _utility_h
#define _utility_h

// read spectrum data from a file

int get_data(double *data, const int bins, const char* filename);


// find regions of interest from the raw data

int find_regions(double *reference, int num_channels, int *regions);


// calculate the area of a peak using Simpson’s rule
// note: returns -1 if the width of the peak is not
// an even number of channels

double calc_peak_area(double *channels, int width);


// Samples can have problems due to interference from
// nearby sources, poor test conditions (temperature, humidity)
// and from deliberate manipulation.

double peak_integrity(double *test, int *region);

// determines if a given test peak matches
// a reference peak subject to the threshold

double match_peak(double *test, double *reference, int *region, double spectrum_energy, double threshold);

#endif