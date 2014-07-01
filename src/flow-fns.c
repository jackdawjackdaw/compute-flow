
#include "flow-fns.h"

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
 * 
 * this doesn't work with two items currently
 */
double subMean(int ns, double prevMean, double xj)
{
  double nmin = ((double)ns) - 1;
  double nr = ((double)ns) / (nmin);
  return(nr * prevMean - (1.0/nmin)*xj);
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

  double ptmin = 0.0;
  double dpt = 0.4;
  
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
