#ifndef __INC_COMPUTE_FLOW__
#define __INC_COMPUTE_FLOW__

#define MAXPTBINS 25


int read_event(FILE* stream, int* ntot, int *ptCount,
                  double ptArray[][MAXPTBINS], double phiArray[][MAXPTBINS],
                  double rapArray[][MAXPTBINS], int chArray[][MAXPTBINS]);


void compute_flow_contrib(int ibin,
                          int *ptCount, int *evCount, double *EventBinFlowContrib,
                         double ptArray[][MAXPTBINS], double phiArray[][MAXPTBINS],
                          double rapArray[][MAXPTBINS], int chArray[][MAXPTBINS]);

void compute_event_plane(int ipart, int ibin,  double *evtPlanes, 
                         int *ptCount,
                         double ptArray[][MAXPTBINS], double phiArray[][MAXPTBINS],
                         double rapArray[][MAXPTBINS], int chArray[][MAXPTBINS]);




double ptmin = 0.0;
double ptmax = 4.0;
double dpt = 0.4;

#endif
