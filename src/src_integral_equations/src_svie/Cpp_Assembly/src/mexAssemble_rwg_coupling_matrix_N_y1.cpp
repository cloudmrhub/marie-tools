#include <iostream>
#include <cmath>
#include <complex>

#include "mex.h"

using std::cout;
using std::endl;
using dcomplex = std::complex <double>;

void Assemble_rwg_coupling_matrix_N_y1(mxComplexDouble * Zbc_Nop, double * Scoord, double * r_p, double * r_n, double * r_2, double * r_3, double * sie_quads, double * vie_quads, double res, const size_t N_vox, const size_t N_rwg, double ko);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{ 
    if(nrhs != 9) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Nine inputs required.");
    }
    
    if(nlhs > 1) {
      mexErrMsgIdAndTxt( "MATLAB:convec:maxlhs",
                         "Too many output arguments.");
    }
    
    double * Scoord    = mxGetDoubles(prhs[0]);
    double * r_p       = mxGetDoubles(prhs[1]);
    double * r_n       = mxGetDoubles(prhs[2]);
    double * r_2       = mxGetDoubles(prhs[3]);
    double * r_3       = mxGetDoubles(prhs[4]);
    double * sie_quads = mxGetDoubles(prhs[5]);
    double * vie_quads = mxGetDoubles(prhs[6]);
    const double res   = mxGetScalar(prhs[7]);
    double ko          = mxGetScalar(prhs[8]);
    
    mwSize N_vox       = mxGetN(prhs[0]);
    mwSize N_rwg       = mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(N_vox, N_rwg, mxCOMPLEX);
    
    mxComplexDouble * Zbc_Gop = mxGetComplexDoubles(plhs[0]);
    
    Assemble_rwg_coupling_matrix_N_y1(Zbc_Gop, Scoord, r_p, r_n, r_2, r_3, sie_quads, vie_quads, res, N_vox, N_rwg, ko);
    
}