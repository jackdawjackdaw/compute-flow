#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// gsl histogram is useful here
#include <gsl/gsl_histogram.h>

#include "compute-v2v3.h"

/**
 * reads input of pt, phi, rap for particles from the awk script pt-phi-grab.
 * computes histograms of v2 and v3 against pt using the event plane method
 * 
 * how is this going to work? we read in all the particles in an event and 
 * store a running average of pt sin(n phi) and pt cos(n phi) for v2, v3
 * 
 * the event plane angles are 
 * 
 * psi_n = (1/n) arctan (  < pt sin (n phi) > / < pt cos (n phi) > ) 
 * 
 * which are easily computed after looping over the particles
 * 
 * finally the averaged flow coeffs are
 * 
 * vn = < cos (n (phi_p - \psi_n ) ) > 
 * 
 * we want to bin these by pt so we compute this for each particle and then fill a histogram
 *
 * so how is everything coming out negative?
 * 
 * 22.10.2013
 * correspondance from HP, 
 * we have to compute the psi_n for each particle, by doing the averaging but 
 * not including the particle we are interested in. 
 *
 * 
 * the code computes a running mean of the event plane angle 
 * for all candidate particles. This mean is then adjusted on a particle by particle basis to remove the 
 * contribution of that particular particle while computing the flow coeffs. 
 * 
 * 
 */




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
  int nevents = 0;

  double ptSinPhi2RM = 0.0; /* running Means for these */
  double ptSinPhi3RM = 0.0; 
  double ptCosPhi2RM = 0.0;
  double ptCosPhi3RM = 0.0;

  /* temp vars */
  double pt;
  double phi;
  double rap;

  double v2var, v3var;
  
  /* histogram stuff */
  double dpt = 0.4;  /* lets make the bin sizes bigger */
  double ptmin = 0.0;
  double ptmax = 4.0;
  int nhistbins = (int)((ptmax - ptmin)/dpt);

  gsl_histogram * v2MeanHist = gsl_histogram_alloc (nhistbins);
  gsl_histogram * v2VarHist = gsl_histogram_alloc (nhistbins);  
  gsl_histogram * v2ptcount = gsl_histogram_alloc (nhistbins);

  gsl_histogram * v3MeanHist = gsl_histogram_alloc (nhistbins);
  gsl_histogram * v3VarHist = gsl_histogram_alloc (nhistbins);  
  gsl_histogram * v3ptcount = gsl_histogram_alloc (nhistbins);

  gsl_histogram_set_ranges_uniform (v2MeanHist, ptmin, ptmax);
  gsl_histogram_set_ranges_uniform (v3MeanHist, ptmin, ptmax);

  gsl_histogram_set_ranges_uniform (v2VarHist, ptmin, ptmax);
  gsl_histogram_set_ranges_uniform (v3VarHist, ptmin, ptmax);
  
  gsl_histogram_set_ranges_uniform (v2ptcount, ptmin, ptmax);
  gsl_histogram_set_ranges_uniform (v3ptcount, ptmin, ptmax);

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
        fprintf(stderr, "#nparts %d nevents %d\n", nparts, nevents); /* this never runs!*/
        for(i = 0; i < nparts; i++){ /* wait, where does this come from? */
          /* we only want to include particles in the FLOW calculation which have rap < 1.0 */
          if(fabs(rapArray[i]) <= rapcutFlow){
            /* just pass in the running means instead of the whole array */
            /* this method is much faster as it avoids a whole extra loop over nparts each time */
            computeFlowRM(2, i, nparts, ptArray[i], phiArray[i], ptSinPhi2RM, ptCosPhi2RM, v2MeanHist, v2VarHist, v2ptcount);
            computeFlowRM(3, i, nparts, ptArray[i], phiArray[i], ptSinPhi3RM, ptCosPhi3RM, v3MeanHist, v3VarHist, v3ptcount);
            //printf("# proc: %d %lf %lf\n", i, ptArray[i], phiArray[i]);
          }
        }

        if(!(nevents % 64))
        /* debug, print the number of events and particles analyzed */
          fprintf(stderr, "# %d %d\n", nevents, nparts);
      }
      // clear the running means
      ptSinPhi2RM = 0.0; 
      ptSinPhi3RM = 0.0; 
      ptCosPhi2RM = 0.0;
      ptCosPhi3RM = 0.0;

      nevents++;
      // have to zero out the number of particles
      nparts = 0;
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
          ptSinPhi3RM = updateMean(nparts, ptSinPhi3RM, pt*sin(3*phi));
          ptCosPhi2RM = updateMean(nparts, ptCosPhi2RM, pt*cos(2*phi));
          ptCosPhi3RM = updateMean(nparts, ptCosPhi3RM, pt*cos(3*phi));

      }
    }

  }

  fprintf(stderr, "#nparts %d nevents %d\n", nparts, nevents); 
  /* push in the final set of particles */
  for(i = 0; i < nparts; i++){
    if(fabs(rapArray[i]) <= rapcutFlow){
      computeFlowRM(2, i, nparts, ptArray[i], phiArray[i], ptSinPhi2RM, ptCosPhi2RM, v2MeanHist, v2VarHist, v2ptcount);
      computeFlowRM(3, i, nparts, ptArray[i], phiArray[i], ptSinPhi3RM, ptCosPhi3RM, v3MeanHist, v3VarHist, v3ptcount);
    }
  }

  printf("# pt dep flow\n");
  printf("# nevents: %d\n", nevents);
  /* output histograms  */
  for(i = 0; i < v2MeanHist->n; i++){
    /* compute the variance directly */
    v2var = (1.0/(v2ptcount->bin[i]))*(v2VarHist->bin[i]) - pow((1.0/(v2ptcount->bin[i]))*v2MeanHist->bin[i], 2.0);
    v3var = (1.0/(v3ptcount->bin[i]))*(v3VarHist->bin[i]) - pow((1.0/(v3ptcount->bin[i]))*v3MeanHist->bin[i], 2.0);

    printf("%d %lf %lf %lf %lf %lf %lf\n", i,
           v2ptcount->bin[i],
           (v2MeanHist->range[i]+v2MeanHist->range[i+1])/2.0,     /* print out the bin centers */
           v2MeanHist->bin[i]/(v2ptcount->bin[i]),
           sqrt(v2var)/sqrt((double)v2ptcount->bin[i]),
           v3MeanHist->bin[i]/v3ptcount->bin[i],
           sqrt(v3var)/sqrt((double)v3ptcount->bin[i])
           );
  }

  gsl_histogram_free (v2MeanHist);
  gsl_histogram_free (v3MeanHist);
  gsl_histogram_free (v2VarHist);
  gsl_histogram_free (v3VarHist);
  gsl_histogram_free (v2ptcount);
  gsl_histogram_free (v3ptcount);
  
  free(ptArray);
  free(phiArray);
  free(rapArray);
  free(buffer);
  
  return EXIT_SUCCESS;
}

/* this increments the mean */
double updateMean(int ns, double prevMean, double newX)
{
  double nplus = (double)ns + 1;
  return( (1.0-(1.0/nplus))*prevMean + (1.0/nplus)*newX);
}

/* given an average: xbar = (1/n) * \sum_{i} Xi  and an Xj that was included in the average,
 *  compute the new average without including Xj.
 * 
 * xbar_{minus} = n/(n-1) xbar - xj / (n-1)
 */
double subMean(int ns, double prevMean, double xj)
{
  double nmin = ((double)ns) - 1;
  double nr = ((double)ns) / (nmin);
  return(nr * prevMean - (1.0/nmin)*xj);
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

