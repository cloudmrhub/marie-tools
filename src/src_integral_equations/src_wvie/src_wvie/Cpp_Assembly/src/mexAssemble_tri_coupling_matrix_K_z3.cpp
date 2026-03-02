#include <iostream>
#include <cmath>
#include <complex>

#include "mex.h"

using std::cout;
using std::endl;
using dcomplex = std::complex <double>;

void Assemble_tri_coupling_matrix_K_z3(mxComplexDouble * Zbc_Nop, double * Scoord, double * r_p, double * r_n, double * r_2, double * wie_quads, double * vie_quads, double res, const size_t N_vox, const size_t N_tri, double ko);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{ 
    if(nrhs != 8) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Eight inputs required.");
    }
    
    if(nlhs > 1) {
      mexErrMsgIdAndTxt( "MATLAB:convec:maxlhs",
                         "Too many output arguments.");
    }
    
    double * Scoord    = mxGetDoubles(prhs[0]);
    double * r_p       = mxGetDoubles(prhs[1]);
    double * r_n       = mxGetDoubles(prhs[2]);
    double * r_2       = mxGetDoubles(prhs[3]);
    double * wie_quads = mxGetDoubles(prhs[4]);
    double * vie_quads = mxGetDoubles(prhs[5]);
    const double res   = mxGetScalar(prhs[6]);
    double ko          = mxGetScalar(prhs[7]);
    
    mwSize N_vox       = mxGetN(prhs[0]);
    mwSize N_tri       = mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(N_vox, N_tri, mxCOMPLEX);
    
    mxComplexDouble * Zbc_Gop = mxGetComplexDoubles(plhs[0]);
    
    Assemble_tri_coupling_matrix_K_z3(Zbc_Gop, Scoord, r_p, r_n, r_2, wie_quads, vie_quads, res, N_vox, N_tri, ko);
    
}