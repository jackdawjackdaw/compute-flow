#ifndef __INC_COMPUTE_V2_NAIIVE__
#define __INC_COMPUTE_V2_NAIIVE__

void computeFlowRM(int n, int id, int nparts,  double pt, double phi,
                   double evtPlaneSinRM, double evtPlaneCosRM, 
                   gsl_histogram* vnMeanHist, gsl_histogram* vnVarHist, gsl_histogram* vnCountHist);

#endif
