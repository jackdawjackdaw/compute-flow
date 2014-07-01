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

void compute_flow_contrib_naiive(int ibin,
                                 int *ptCount, int *evCount, double *EventBinFlowContrib,
                                 double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                                 double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);



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
  int nevents = 0;

  double meanV2, varV2;
  double nev;
  
  double binCent = 0.0;
  
  /* histogram stuff */
  double dpt = 0.4;  /* lets make the bin sizes bigger */
  double ptmin = 0.0;
  double ptmax = 4.0;

  double v2MeanHist[MAXPTBINS];
  double v2SqHist[MAXPTBINS];

  for(i = 0; i < MAXPTBINS; i++){
    v2MeanHist[i] = 0.0;
    v2SqHist[i] = 0.0;
    evCount[i] = 0.0;
  }
  
  
  
  char * buffer = NULL;
  size_t linecap = 0;
  size_t linelen;
  size_t retval;

  double EventBinFlow[2] = {0.0, 0.0};
  
  /* get the first line, which is always an evt line, so we can ignore */
  getline(&buffer, &linecap, stdin);

  do {
    reset_arrays(ptCount, ptArray, phiArray, rapArray, chArray);
    
    retval = read_event(stdin, &nparts, ptCount, ptArray, phiArray, rapArray, chArray);
    printf("# %d %d\n\r", nevents, nparts);
    nevents++;

    for(i = 0; i < MAXPTBINS; i++){
      compute_flow_contrib_naiive(i, ptCount, evCount, EventBinFlow,
                                  ptArray, phiArray, rapArray, chArray);

      v2MeanHist[i] += EventBinFlow[0];
      v2SqHist[i]  += EventBinFlow[1];
    }
    
  } while (retval != EOF); /* we finished the stream */

    for(i = 0; i < MAXPTBINS; i++){
    if(evCount[i] > 0){
      v2MeanHist[i] /= (double)evCount[i];
      v2SqHist[i]  /= (double)evCount[i];
    }
  }

  /* finally we can print everything out */
  for( i = 0 ; i < MAXPTBINS; i++){
    binCent = (double)i*dpt + ptmin + (dpt)/2;
    nev = (double)evCount[i];
    if(evCount[i] > 1 && binCent <= (ptmax+0.5)){
      meanV2 = v2MeanHist[i];

      varV2 = ((v2SqHist[i])-meanV2*meanV2);

      printf("%lf %d %lf %lf\n",
             binCent, evCount[i],
             meanV2, sqrt(varV2)/(nev-1));
    }
  }

  
  free(buffer);
  
  return EXIT_SUCCESS;
}




void compute_flow_contrib_naiive(int ibin,
                          int *ptCount, int *evCount, double *EventBinFlowContrib,
                          double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                          double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS])
{

  int i,nch=0;
  double rapCutFlow = 1.0;
  double v2Sum = 0.0;
  double px, py;

  if(ptCount[ibin] > 0){
    evCount[ibin]++; /* this event has contributed to this bin*/
    for(i = 0; i < ptCount[ibin]; i++){
      if(abs(chArray[i][ibin])>0 && abs(rapArray[i][ibin]) <= rapCutFlow){
        nch++;
        /* do we really need to do this each time? */

        px = ptArray[i][ibin]*cos(phiArray[i][ibin]);
        py = ptArray[i][ibin]*sin(phiArray[i][ibin]);

        /* do we want to average the whole thing or the top and bottom separately?
         * start by looking at the whole thing */
        v2Sum += (px*px - py*py ) / (px*px + py*py);
        //v2Sum += cos(2*phiArray[i][ibin]);
      }
    }
  }

  
  if(nch > 0){
    EventBinFlowContrib[0] = v2Sum/(double)nch;
    EventBinFlowContrib[1] = pow(v2Sum/(double)nch, 2.0);
  } else {
    EventBinFlowContrib[0] = 0.0;
    EventBinFlowContrib[1] = 0.0;
  }
  
}
  
  
