#ifndef __INC_COMPUTE_FLOW__
#define __INC_COMPUTE_FLOW__
#define MAXPARTS 10000
#define MAXPTBINS 25

void reset_arrays(int *ptCount,
                  double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                  double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);


int read_event(FILE* stream, int* ntot, int *ptCount,
                  double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                  double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);


void compute_flow_contrib(int ibin,
                          int *ptCount, int *evCount, double *EventBinFlowContrib,
                         double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                          double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);

void compute_event_plane(int ipart, int ibin,  double *evtPlanes, 
                         int *ptCount,
                         double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                         double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);




double ptmin = 0.0;
double ptmax = 4.0;
double dpt = 0.4;

#endif
