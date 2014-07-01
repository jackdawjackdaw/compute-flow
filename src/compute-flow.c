#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compute-flow.h"

/**
 * ccs, cec24@phy.duke.edu, 25.06.2014
 *
 * a retry at computing the flow coeffs, start by explicitly matching HPs calculation
 */


int main (int argc, char* argv[]){
  int i;
  /* double rapcutPlane = 2.0; /\* just for testing, change this *\/ */
  /* double rapcutFlow = 1.0; */
  double binCent = 0.0;

  /* contains pt counts */
  int ptCount[MAXPTBINS];
  int evCount[MAXPTBINS];
  /* holds the pt etc of the particles in the event, indexed by the count (in that bin) and the bin */
  double ptArray[MAXPARTS][MAXPTBINS];
  double phiArray[MAXPARTS][MAXPTBINS];
  double rapArray[MAXPARTS][MAXPTBINS];
  int chArray[MAXPARTS][MAXPTBINS];


  double v2MeanHist[MAXPTBINS];
  double v2SqHist[MAXPTBINS];
  double v3MeanHist[MAXPTBINS];
  double v3SqHist[MAXPTBINS];

  double meanV2, meanV3, varV2, varV3;
  double nev;
  
  int nparts = 0;
  int nevents = 0;
  int retval = 0;

  char * buffer = NULL;
  size_t linecap = 0;
  size_t linelen;

  double EventBinFlow[4];
  
  for(i = 0; i < MAXPTBINS; i++){
    v2MeanHist[i] = 0.0;
    v2SqHist[i] = 0.0;
    v3MeanHist[i] = 0.0;
    v3SqHist[i] = 0.0;
    evCount[i] = 0;
  }


  /* get the first line, which is always an evt line, so we can ignore */
  getline(&buffer, &linecap, stdin);
  
  do {
    reset_arrays(ptCount, ptArray, phiArray, rapArray, chArray);
    
    retval = read_event(stdin, &nparts, ptCount, ptArray, phiArray, rapArray, chArray);
    nevents++;
    //printf("# %d %d\n\r", nevents, nparts);

    if(nevents % 64 == 0)
      fprintf(stderr, "# %d %d\n\r", nevents, nparts);
    
    /* process the event into the bins */
    for(i = 0; i < MAXPTBINS; i++){

      EventBinFlow[0] = 0.0;
      EventBinFlow[1] = 0.0;
      EventBinFlow[2] = 0.0;
      EventBinFlow[3] = 0.0;
      
      compute_flow_contrib(i, ptCount, evCount, EventBinFlow,
                           ptArray, phiArray, rapArray, chArray);

      v2MeanHist[i] += EventBinFlow[0];
      v2SqHist[i]  += EventBinFlow[1];
      v3MeanHist[i] += EventBinFlow[2];
      v3SqHist[i]  += EventBinFlow[3];

      /* printf("%d %lf %lf %lf %lf\n", i, */
      /*        EventBinFlow[0], EventBinFlow[1], EventBinFlow[2], EventBinFlow[3]); */
    }
      
    
  } while (retval != EOF); /* we finished the stream */

  for(i = 0; i < MAXPTBINS; i++){
    if(evCount[i] > 0){
      v2MeanHist[i] /= (double)evCount[i];
      v2SqHist[i]  /= (double)evCount[i];
      v3MeanHist[i] /= (double)evCount[i];
      v3SqHist[i]  /= (double)evCount[i];
    }
  }

  /* finally we can print everything out */
  for( i = 0 ; i < MAXPTBINS; i++){
    binCent = (double)i*dpt + ptmin + (dpt)/2;
    nev = (double)evCount[i];
    if(evCount[i] > 1){
      meanV2 = v2MeanHist[i];
      meanV3 = v3MeanHist[i];

      varV2 = (v2SqHist[i])-meanV2*meanV2;
      varV3 = (v3SqHist[i])-meanV3*meanV3;

      printf("%lf %d %lf %lf %lf %lf\n",
             binCent, evCount[i],
             meanV2, sqrt(varV2)/(nev-1),
             meanV3, sqrt(varV3)/(nev-1));
    }
  }
      
}

void compute_flow_contrib(int ibin,
                          int *ptCount, int *evCount, double *EventBinFlowContrib,
                         double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                         double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS])
  
{
  int nch = 0;
  int i;
  double rapCutFlow = 1.0;
  double evtPlanes[2] = {0.0, 0.0};

  double v2Sum = 0.0;
  double v3Sum = 0.0;

  if(ptCount[ibin] > 0){
    evCount[ibin]++; /* this event has contributed to this bin*/
    for(i = 0; i < ptCount[ibin]; i++){
      if(abs(chArray[i][ibin])>0 && abs(rapArray[i][ibin]) <= rapCutFlow){
        nch++;
        compute_event_plane(i, ibin, evtPlanes,
                            ptCount, ptArray, phiArray, rapArray, chArray);

        v2Sum += cos(2*phiArray[i][ibin] - evtPlanes[0]);
        v3Sum += cos(3*phiArray[i][ibin] - evtPlanes[1]);
      }
    }
  }

  if(nch > 0){
    EventBinFlowContrib[0] = v2Sum/(double)nch;
    EventBinFlowContrib[1] = pow(v2Sum/(double)nch, 2.0);
    EventBinFlowContrib[2] = v3Sum/(double)nch;
    EventBinFlowContrib[3] = pow(v3Sum/(double)nch, 2.0);
  } else {
    EventBinFlowContrib[0] = 0.0;
    EventBinFlowContrib[1] = 0.0;
    EventBinFlowContrib[2] = 0.0;
    EventBinFlowContrib[3] = 0.0;
  }
}

/* we loop over all the particles not in this pt bin */
void compute_event_plane(int ipart, int ibin,  double *evtPlanes, 
                         int *ptCount,
                         double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                         double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS])
{
  int i, j;
  int nloop = 0;
  double v2SinSum = 0.0;
  double v2CosSum = 0.0;
  double v3SinSum = 0.0;
  double v3CosSum = 0.0;

  double rapCutEvtPlane = 2.0;

  for(i = 0; i < MAXPTBINS; i++){
    if(ptCount[i] > 0){

      for(j = 0; j < ptCount[i]; j++){
        if(abs(rapArray[j][i]) <= rapCutEvtPlane){
          if(i==ibin && j == ipart){
            continue;
          }
          nloop++;

          v2SinSum += ptArray[j][i]*sin(2*phiArray[j][i]);
          v2CosSum += ptArray[j][i]*cos(2*phiArray[j][i]);

          v3SinSum += ptArray[j][i]*sin(3*phiArray[j][i]);
          v3CosSum += ptArray[j][i]*cos(3*phiArray[j][i]);
        }
      }
    }
  }

  v2SinSum = v2SinSum /(double)nloop;
  v2CosSum = v2CosSum /(double)nloop;
  v3SinSum = v3SinSum /(double)nloop;
  v3CosSum = v3CosSum /(double)nloop;

  evtPlanes[0] = atan2(v2SinSum, v2CosSum)/2.0;
  evtPlanes[1] = atan2(v3SinSum, v3CosSum)/3.0;
  
}

void reset_arrays(int *ptCount,
                  double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                  double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS])
{
  int i,j;
  for(i = 0; i < MAXPTBINS; i++) ptCount[i] = 0;

  for(i = 0; i < MAXPARTS; i++)
    for(j = 0; j < MAXPTBINS; j++){
      ptArray[i][j]  = 0.0;
      phiArray[i][j] = 0.0;
      chArray[i][j]  = 0;
    }
  
}

int read_event(FILE* stream, int* ntot, int *ptCount,
                  double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                  double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS])
{
  char * buffer = NULL;
  size_t linecap = 0;
  size_t linelen;

  int charge;
  double pt, phi, rap;
  int nparts = 0;
  int ipt, pcount;
  
  while((linelen = getline(&buffer, &linecap, stream)) != EOF){
    
    if(strncmp(buffer, "#", 1) == 0){ /* break on #event lines */
      /* store the total number read */
      *ntot = nparts;
      return 0;
    }

    sscanf(buffer, "%d %lf %lf %lf",&charge,  &pt, &phi, &rap);
    nparts++;

    //printf("%d %lf %lf %lf\n", charge, pt, phi, rap);
    
    /* compute which pt bin this particle goes into */
    /* this rounding is perhaps a source of pt binning confusion */
    //ipt = (lround)((pt-ptmin)/dpt);
    ipt = (int)((pt-ptmin)/dpt);

    if(ipt < MAXPTBINS){
      ptCount[ipt]++;
      pcount = ptCount[ipt];
      ptArray[pcount][ipt]  = pt;
      phiArray[pcount][ipt] = phi;
      rapArray[pcount][ipt] = rap;
      chArray[pcount][ipt]  = charge;
    }

  }

  
  /* store the total number read */
  *ntot = nparts;

  
  if(linelen == EOF)
    return EOF;

  /* shouldn't ever reach here */
}
