#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "flow-fns.h"
#include "compute-flow.h"


void compute_evt_plane_avg(double* evtPlaneMeans,
                           int* nEvtPlane,
                           int* ptCount,
                           double ptArray[][MAXPTBINS], double phiArray[][MAXPTBINS],
                           double rapArray[][MAXPTBINS], int chArray[][MAXPTBINS]);

void compute_flow_contrib_rm(int ibin,
                             int *ptCount, int *evCount, double *EventBinFlowContrib,
                             double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                             double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS]);
  


/**
 * ccs, cec24@phy.duke.edu, 25.06.2014
 *
 * a retry at computing the flow coeffs, start by explicitly matching HPs calculation
 *
 * doesn't perfectly agree with the output from HP's calc
 * also it's rather slow
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

  /* set the bin width in gev */
  double dpt = 0.4;
  double ptmin = 0.0;
  double ptmax = 4.0;
  
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
    
    retval = read_event(stdin, dpt,&nparts, ptCount, ptArray, phiArray, rapArray, chArray);
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
      
      compute_flow_contrib_rm(i, ptCount, evCount, EventBinFlow,
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
    if(evCount[i] > 1 && binCent <= (ptmax+0.5)){
      meanV2 = v2MeanHist[i];
      meanV3 = v3MeanHist[i];

      varV2 = ((v2SqHist[i])-meanV2*meanV2);
      varV3 = (v3SqHist[i])-meanV3*meanV3;

      printf("%lf %d %lf %lf %lf %lf\n",
             binCent, evCount[i],
             meanV2, sqrt(varV2)/(nev-1),
             meanV3, sqrt(varV3)/(nev-1));
    }
  }

  free(buffer);
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
        /* do we really need to do this each time? */
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


/* agrees with old method
 * old method runs in 37seconds on the test set
 * this method runs in 0.29seconds on the same test set
 */
void compute_flow_contrib_rm(int ibin,
                          int *ptCount, int *evCount, double *EventBinFlowContrib,
                         double ptArray[MAXPARTS][MAXPTBINS], double phiArray[MAXPARTS][MAXPTBINS],
                         double rapArray[MAXPARTS][MAXPTBINS], int chArray[MAXPARTS][MAXPTBINS])
  
{
  int nch = 0;
  int i;
  int nPartsEvtPlane = 0;
  double rapCutFlow = 1.0;
  double evtPlanesAvg[4] = {0.0, 0.0, 0.0, 0.0};
  double evtPlanes[2] = {0.0, 0.0};

  double meanCosVn = 0.0;
  double meanSinVn = 0.0;

  double cosVn = 0.0;
  double sinVn = 0.0;
  
  double v2Sum = 0.0;
  double v3Sum = 0.0;

  if(ptCount[ibin] > 0){
    evCount[ibin]++; /* this event has contributed to this bin*/
    
    /* first compute the event plane */
    compute_evt_plane_avg(evtPlanesAvg, &nPartsEvtPlane, ptCount, ptArray, phiArray, rapArray, chArray);
    
    for(i = 0; i < ptCount[ibin]; i++){
      if(abs(chArray[i][ibin])>0 && abs(rapArray[i][ibin]) <= rapCutFlow){
        nch++;

        /* now we need to restrict our evtPlanes to the ones without this
         * particular particle in it */
        cosVn = ptArray[i][ibin] * cos(2.0*phiArray[i][ibin]);
        sinVn = ptArray[i][ibin] * sin(2.0*phiArray[i][ibin]);
        
        meanCosVn = subMean(nPartsEvtPlane, evtPlanesAvg[0], cosVn );
        meanSinVn = subMean(nPartsEvtPlane, evtPlanesAvg[1], sinVn);        
        evtPlanes[0] = atan2(meanSinVn, meanCosVn)/2.0;

        /* again for v3 */
        cosVn = ptArray[i][ibin] * cos(3.0*phiArray[i][ibin]);
        sinVn = ptArray[i][ibin] * sin(3.0*phiArray[i][ibin]);

        meanCosVn = subMean(nPartsEvtPlane, evtPlanesAvg[2], cosVn);
        meanSinVn = subMean(nPartsEvtPlane, evtPlanesAvg[3], sinVn);        
        evtPlanes[1] = atan2(meanSinVn, meanCosVn)/3.0;

        /* now compute the actual contribution of this particle */
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

/** 
 * compute the event plane (cos and sin) averaged over all the bins in the event 
 */
void compute_evt_plane_avg(double* evtPlaneMeans, int* nEvtPlane,
                           int* ptCount,
                           double ptArray[][MAXPTBINS], double phiArray[][MAXPTBINS],
                           double rapArray[][MAXPTBINS], int chArray[][MAXPTBINS])
{

  int i, j;
  int nloop = 0;

  double meanSinV2 = 0.0;
  double meanCosV2 = 0.0;
  double meanSinV3 = 0.0;
  double meanCosV3 = 0.0;

  
  double rapCutEvtPlane = 2.0;
  
  for(i = 0; i < MAXPTBINS; i++){
    for(j = 0; j < ptCount[i]; j++){
      if(abs(rapArray[j][i]) <= rapCutEvtPlane){
        nloop++;
        meanSinV2 += ptArray[j][i]*sin(2.0*phiArray[j][i]);
        meanCosV2 += ptArray[j][i]*cos(2.0*phiArray[j][i]);

        meanSinV3 += ptArray[j][i]*sin(3.0*phiArray[j][i]);
        meanCosV3 += ptArray[j][i]*cos(3.0*phiArray[j][i]);

      }
    }
  }

  meanCosV2 = meanCosV2 / (double)nloop;
  meanCosV3 = meanCosV3 / (double)nloop;

  meanSinV2 = meanSinV2 / (double)nloop;
  meanSinV3 = meanSinV3 / (double)nloop;  

  evtPlaneMeans[0] = meanCosV2;
  evtPlaneMeans[1] = meanSinV2;
  evtPlaneMeans[2] = meanCosV3;
  evtPlaneMeans[3] = meanSinV3;


  /* store the number of particles which contributed */
  *nEvtPlane = nloop;
}




/* we loop over all the particles not in this pt bin
 * this is decidedly slow, we don't need to recompute the average each time
 */
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

