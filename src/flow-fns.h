#ifndef __flow_fns__
#define __flow_fns__

#define MAXPARTS 10000
#define MAXPTBINS 25

#include <stdio.h>
#include <string.h>

double updateMean(int ns, double prevMean, double newX);
double subMean(int ns, double prevMean, double xj);

void reset_arrays(int *ptCount,
                  double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                  double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);

int read_event(FILE* stream, double dpt, int* ntot, int *ptCount, 
               double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
               double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);


#endif
