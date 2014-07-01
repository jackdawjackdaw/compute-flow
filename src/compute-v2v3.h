#ifndef __INC_COMPUTE_V2V3__
#define __INC_COMPUTE_V2V3__

void computeFlow(int n, int id, int nparts,  double pt, double phi,
                 double* evtPlaneSinArray, double* evtPlaneCosArray, 
                 gsl_histogram* vnMeanHist, gsl_histogram* vnVarHist, gsl_histogram* vnCountHist);

void computeFlowRM(int n, int id, int nparts,  double pt, double phi,
                   double evtPlaneSinRM, double evtPlaneCosRM, 
                   gsl_histogram* vnMeanHist, gsl_histogram* vnVarHist, gsl_histogram* vnCountHist);

#endif
