
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
