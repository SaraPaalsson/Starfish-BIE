#include "mex.h"
#include <math.h>
#include <string.h>
double times_two(double W0[]){
    double A;
    return A = 2*W0[0];
}
void mexFunction(int nlhs, mxArray * plhs[], int nrhs,
                 const mxArray * prhs[]) {
    double A, *W0;
    double *tmp;
    mwIndex idx = 1;
    int ifield, nfields;
    const char **fname;       /* pointers to field names */
      nfields = mxGetNumberOfFields(prhs[0]);
      fname = mxCalloc(nfields, sizeof(*fname));
      for (ifield=0; ifield< nfields; ifield++){
          fname[ifield] = mxGetFieldNameByNumber(prhs[0],ifield);
          if (strcmp(fname[ifield],"W0") == 0){
              tmp = mxGetPr(mxGetFieldByNumber(prhs[0],idx,ifield));
              W0 = tmp;
          }
      }
      A = times_two(W0);
      mexPrintf("A = %%f",A);
  }