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
  double rapcutPlane = 2.0; /* just for testing, change this */
  double rapcutFlow = 1.0;

  double* ptArray;
  double* phiArray;
  double* rapArray;

  int nparts = 0;
  int nch = 0;
  int nevents = 0;

  double ptSinPhi2RM = 0.0; /* running Means for these */
  double ptCosPhi2RM = 0.0;

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

  gsl_histogram *v2EvtHist = gsl_histogram_alloc (nhistbins);
  gsl_histogram *v2EvtVarHist = gsl_histogram_alloc (nhistbins);
  gsl_histogram *v2EvtCount = gsl_histogram_alloc (nhistbins);

  gsl_histogram_set_ranges_uniform (v2EvtHist, ptmin, ptmax);
  gsl_histogram_set_ranges_uniform (v2EvtVarHist, ptmin, ptmax);
  gsl_histogram_set_ranges_uniform (v2EvtCount, ptmin, ptmax);
  
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

  ptArray = malloc(sizeof(double) * MAXPARTS); 
  phiArray = malloc(sizeof(double) * MAXPARTS);
  rapArray = malloc(sizeof(double) * MAXPARTS);

  while( (linelen = getline(&buffer, &linecap, stdin)) != -1) { 
    if(strncmp(buffer, "#", 1) == 0){ /* break on #event lines */
      if(nevents > 0){
        /* compute the flow contribution of the particles for this event*/
        //fprintf(stderr, "#nparts %d nevents %d\n", nparts, nevents); 
        for(i = 0; i < nparts; i++){ 
          /* we only want to include particles in the FLOW calculation which have rap < 1.0 */
          if(fabs(rapArray[i]) <= rapcutFlow){
            /* just pass in the running means instead of the whole array */
            /* this method is much faster as it avoids a whole extra loop over nparts each time */
            nch++;
            computeFlowRM(2, i, nparts, ptArray[i], phiArray[i], ptSinPhi2RM, ptCosPhi2RM, v2EvtHist, v2EvtVarHist, v2EvtCount);
            //printf("# proc: %d %lf %lf\n", i, ptArray[i], phiArray[i]);
          }
        }

        if(!(nevents % 64))
        /* debug, print the number of events and particles analyzed */
          fprintf(stderr, "# %d %d\n", nevents, nparts);
      }
      // clear the running means
      ptSinPhi2RM = 0.0; 
      ptCosPhi2RM = 0.0;

      nevents++;

      /* now increment the actual means*/
      for(i = 0; i < nhistbins; i++){
          v2Mean[i] += (v2EvtHist->bin[i]);///v2EvtCount->bin[i]);
          v2Var[i] += pow(v2EvtVarHist->bin[i], 2.0);
          v2Count[i]++;
      }
      /* now we zero the event histograms */
      //gsl_histogram_reset(v2EvtCount);
      gsl_histogram_reset(v2EvtHist);
      gsl_histogram_reset(v2EvtVarHist);
      
      // have to zero out the number of particles
      nparts = 0;
      nch = 0;

      
    } else {  /* actually read the particle lines from an event */
      // split the buffer into pt, phi and rap
      // store pt and phi and event plane contributions for each particls
      sscanf(buffer, "%lf %lf %lf", &pt, &phi, &rap);
      if(fabs(rap) <= rapcutPlane && phi != 0.0 && rap != 0.0 ){
          ptArray[nparts] = pt;
          phiArray[nparts] = phi;
          rapArray[nparts] = rap;
          //printf("# init: %d %lf %lf\n", nparts, ptArray[nparts], phiArray[nparts]);
          nparts++;

          // increment the running means
          ptSinPhi2RM = updateMean(nparts, ptSinPhi2RM, pt*sin(2*phi));
          ptCosPhi2RM = updateMean(nparts, ptCosPhi2RM, pt*cos(2*phi));

      }
    }

  }

  //fprintf(stderr, "#nparts %d nevents %d\n", nparts, nevents); 
  /* push in the final set of particles */
  /* lets just ignore this for a moment */
  for(i = 0; i < nparts; i++){
    if(fabs(rapArray[i]) <= rapcutFlow){
      computeFlowRM(2, i, nparts, ptArray[i], phiArray[i], ptSinPhi2RM, ptCosPhi2RM, v2EvtHist, v2EvtVarHist, v2EvtCount);
    }
  }

  for(i = 0; i < nhistbins; i++){
      v2Mean[i] += (v2EvtHist->bin[i]);
      v2Var[i] += pow(v2EvtVarHist->bin[i], 2.0);
      v2Count[i]++;
  }

  for(i = 0; i < nhistbins; i++){
    if(v2EvtCount->bin[i] > 0){
      v2Mean[i] = v2Mean[i] / (double)v2EvtCount->bin[i];
      v2Var[i] = v2Var[i] / pow((double)v2EvtCount->bin[i], 2.0);
    }
    //printf("bar %d %lf %lf\n", i, v2EvtCount->bin[i], v2Mean[i]);    
  }

  printf("# pt dep flow\n");
  printf("# nevents: %d\n", nevents);
  /* output histograms  */
  for(i = 0; i < nhistbins; i++){
    /* compute the variance directly */
    if(v2Count[i] > 0){
      /* v2mean = (1.0/((double)v2Count[i]))*v2Mean[i]; */
      /* v2var = (1.0/((double)v2Count[i]))*(v2Var[i]) - pow(v2mean, 2.0); */

      v2mean = v2Mean[i];
      v2var = (v2Var[i]) - pow(v2mean, 2.0);

    } else {
      v2mean = 0.0;
      v2var = 0.0;
    }

    printf("%d %d %lf %lf %lf\n", i,
           v2Count[i],
           (v2EvtHist->range[i]+v2EvtHist->range[i+1])/2.0,     /* print out the bin centers */
           v2mean,
           sqrt(v2var/((double)v2Count[i]-1))
           );
  }

  free(v2Mean);
  free(v2Var);
  free(v2Count);

  gsl_histogram_free (v2EvtHist);
  gsl_histogram_free (v2EvtVarHist);
  gsl_histogram_free (v2EvtCount);

  
  free(ptArray);
  free(phiArray);
  free(rapArray);
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

