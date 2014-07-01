#ifndef __flow_fns__
#define __flow_fns__
#define MAXPARTS 10000
#define MAXPTBINS 25

#include <stdio.h>
#include <string.h>

double updateMean(int ns, double prevMean, double newX);
double subMean(int ns, double prevMean, double xj);

int read_event(FILE* stream, int* ntot, int *ptCount,
                  double ptArray[][MAXPTBINS], double phiArray[][MAXPTBINS],
                 double rapArray[][MAXPTBINS], int chArray[][MAXPTBINS]);


void reset_arrays(int *ptCount,
                  double ptArray[][MAXPTBINS], double phiArray[][MAXPTBINS],
                    double rapArray[][MAXPTBINS], int chArray[][MAXPTBINS]);

#endif
