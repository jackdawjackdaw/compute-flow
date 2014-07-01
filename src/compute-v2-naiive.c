#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// gsl histogram is useful here
#include <gsl/gsl_histogram.h>

/**
 * reads input of pt, phi, rap for particles from the awk script pt-phi-grab.
 * computes histogram of v2 against pt using the naiive definition, well it should do this, it doesnt' yet
 * 
 * v2 = < px**2 - py**2 > / < px**2 + py**2>
 * 
 * \todo this is totally unfinished
 */

#include "compute-v2-naiive.h"

#include "flow-fns.h"


/**
 * should read command line options to set the pt bin size and the dy window 
 * for now we'll fix this to being < 1
 */
int main (int argc, char* argv[]){
  int i;
  double rapcutFlow = 1.0;

  int ptCount[MAXPTBINS];
  int evCount[MAXPTBINS];
  double ptArray[MAXPARTS][MAXPTBINS];
  double phiArray[MAXPARTS][MAXPTBINS];
  double rapArray[MAXPARTS][MAXPTBINS];
  int chArray[MAXPARTS][MAXPTBINS];

  int nparts = 0;
  int nch = 0;
  int nevents = 0;

  /* temp vars */
  double pt;
  double phi;
  double rap;

  double v2var, v3var;
  double v2mean;
  
  /* histogram stuff */
  double dpt = 0.4;  /* lets make the bin sizes bigger */
  double ptmin = 0.0;
  double ptmax = 4.0;
  int nhistbins = (int)((ptmax - ptmin)/dpt);

  double *v2Mean = malloc(sizeof(double)*nhistbins);
  double *v2Var = malloc(sizeof(double)*nhistbins);
  int *v2Count = malloc(sizeof(int)*nhistbins);  
  
  for(i = 0; i < nhistbins; i++){
    v2Mean[i] = 0.0;
    v2Var[i] = 0.0;
    v2Count[i] = 0;
  }
  
  char * buffer = NULL;
  size_t linecap = 0;
  size_t linelen;
  size_t retval;

  /* get the first line, which is always an evt line, so we can ignore */
  getline(&buffer, &linecap, stdin);

  do {
    reset_arrays(ptCount, ptArray, phiArray, rapArray, chArray);
    
    retval = read_event(stdin, &nparts, ptCount, ptArray, phiArray, rapArray, chArray);
    printf("# %d %d\n\r", nevents, nparts);
    nevents++;
  } while (retval != EOF); /* we finished the stream */


  free(v2Mean);
  free(v2Var);
  free(v2Count);
  free(buffer);
  
  return EXIT_SUCCESS;
}



void computeFlow(int n, int id, int nparts,  double pt, double phi, double* evtPlaneSinArray, double* evtPlaneCosArray, 
                 gsl_histogram* vnMeanHist, gsl_histogram* vnVarHist, gsl_histogram* vnCountHist)
{
  /* compute the event plane angle for this particle 
  *  importantly this includes the contribution from all the other particles but this one
  * start with a naiive algorithm for this
  * 
  * well this is silly, if we have the average, we can just correct it by removing the contribution of the particle we care about
  * right?
  */
  int i;
  double psi = 0.0;
  double meanSin = 0.0, meanCos = 0.0;
  double flowCoeff = 0.0;
  for(i = 0; i < nparts; i ++){
    if( i != id){
      meanSin += evtPlaneSinArray[i];
      meanCos += evtPlaneCosArray[i];
    }
  }
  meanSin = (1.0/(double)nparts)*meanSin;
  meanCos = (1.0/(double)nparts)*meanCos;

  /* the event plane angle for *THIS* particle */
  psi = (1.0/(double)n) * atan2(meanSin, meanCos);
    
  /* now accumulate the histograms with the right things */
  flowCoeff = cos( ((double)n) * (phi - psi));

  /* if(id < 10 && n == 2){ */
  /*   printf("%d %lf %lf %lf %lf\n", id, */
  /*          meanSin, meanCos, psi, flowCoeff); */
  /* } */

  
  gsl_histogram_accumulate(vnMeanHist, pt, flowCoeff);
  gsl_histogram_accumulate(vnVarHist, pt, flowCoeff*flowCoeff);  
  gsl_histogram_increment(vnCountHist, pt);
}


void computeFlowRM(int n, int id, int nparts,  double pt, double phi, double evtPlaneSinRM, double evtPlaneCosRM, 
                 gsl_histogram* vnMeanHist, gsl_histogram* vnVarHist, gsl_histogram* vnCountHist)
{
  /* compute the event plane angle for this particle 
  *  importantly this includes the contribution from all the other particles but this one
  * start with a naiive algorithm for this
  * 
  * well this is silly, if we have the average, we can just correct it by removing the contribution of the particle we care about
  * right?
  */
  double psi = 0.0;
  double meanSin = 0.0, meanCos = 0.0;
  double flowCoeff = 0.0;

  meanSin = subMean(nparts, evtPlaneSinRM, pt*sin((double)n*phi));
  meanCos = subMean(nparts, evtPlaneCosRM, pt*cos((double)n*phi));  
  
  /* the event plane angle for *THIS* particle */
  psi = (1/(double)n) * atan2(meanSin, meanCos);
    
  /* now accumulate the histograms with the right things */
  flowCoeff = cos( ((double)n) * (phi - psi));

  /* if(id < 10 && n == 2){ */
  /*   printf("%d %lf %lf %lf %lf\n", id, */
  /*          meanSin, meanCos, psi, flowCoeff); */
  /* } */

  gsl_histogram_accumulate(vnMeanHist, pt, flowCoeff);
  gsl_histogram_accumulate(vnVarHist, pt, flowCoeff*flowCoeff);  
  gsl_histogram_increment(vnCountHist, pt);
}

