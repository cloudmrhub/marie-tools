#include <iostream>
#include <iomanip>
#include <math.h>
#include <complex>
#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using std::cout;
using std::endl;
using dcomplex = std::complex <double>;
using namespace std;
  
void get_source_coords_mat(double * r_v, double * r_2, double * wie_quads, double * z_wie, double * t, const size_t Np_wie, const size_t N_tri);
double vec_norm_l2(double * r);
double* compute_segment_lengths(double * r_p, double * r_2, const size_t N_tri);

///////////////////////////////////////////////////////////////////////////
void  Assemble_tri_coupling_matrix_N_x3(mxComplexDouble * Zbc_Nop, double * Scoord, double * r_p, double * r_n, double * r_2, double * wie_quads, double * vie_quads, double res, const size_t N_vox, const size_t N_tri, double ko)
{
    const size_t N2 = 2;
    const size_t N3 = 3;
    const size_t Np_wie = wie_quads[0];  
    const size_t Np_vie = vie_quads[0];

    // get scaling constants 
    double mu_0 = 4 * M_PI * 1e-7;
    double c_0  = 299792458;
    double e_0  = 1.0 / c_0 / c_0 / mu_0; 
    dcomplex ce = dcomplex(0.0,1.0) * ko * c_0 * e_0;
    dcomplex em_scaling = - 1.0 / ce / 4.0 / M_PI;
    dcomplex m_scaling = - 1.0 / 4.0 / M_PI;
            
    // allocate memory for quadrature points for positive and negative lines 
    double * r_src_p = new double[N3 * Np_wie * N_tri];
    double * r_src_n = new double[N3 * Np_wie * N_tri];
    double * f_wie_p = new double[N3 * Np_wie * N_tri];
    double * f_wie_n = new double[N3 * Np_wie * N_tri];
    
    // get coordinates of quadrature points on the segments of a given line pair
    get_source_coords_mat(r_2, r_p, wie_quads, r_src_p, f_wie_p, Np_wie, N_tri);
    get_source_coords_mat(r_n, r_2, wie_quads, r_src_n, f_wie_n, Np_wie, N_tri);   


    // Compute segment lengths for positive and negative segments
    double * L_p = compute_segment_lengths(r_p, r_2, N_tri);
    double * L_n = compute_segment_lengths(r_n, r_2, N_tri);
    
    // loop over the voxels (this part should run in parallel)
#ifdef _OPENMP
   #pragma omp parallel for
#endif
    for (size_t i_vox = 0; i_vox < N_vox; ++i_vox) {
        
        // allocate space for observer coordinates
        double r_obs_voxel[N3];

        // get center of a current voxel
        double r_vox_center[N3] = {Scoord[0 + i_vox * N3],
                                   Scoord[1 + i_vox * N3],
                                   Scoord[2 + i_vox * N3]};
        
        // loop other quadrature points along x axis
        for (size_t i_vie = 0; i_vie < Np_vie; ++i_vie) {
            
            // define x-component of the observation point 
            r_obs_voxel[0] = r_vox_center[0] + res / 2.0 * vie_quads[i_vie + Np_vie + 1];
            
            // loop other quadrature points along y axis
            for (size_t j_vie = 0; j_vie < Np_vie; ++j_vie) {
                
                // define y-component of the observation point 
                r_obs_voxel[1] = r_vox_center[1] + res / 2.0 * vie_quads[j_vie + Np_vie + 1];
                
                for (size_t k_vie = 0; k_vie < Np_vie; ++k_vie) {
                    
                    // sdefine z-component of the observation point
                    r_obs_voxel[2] = r_vox_center[2] + res / 2.0 * vie_quads[k_vie + Np_vie + 1];
                
                    // allocate memory for the vector of distances 
                    double R_vec [N3];
                    
                    for (size_t i_tri = 0; i_tri < N_tri; ++i_tri) {
                        
                        size_t tri_shift_in  = i_tri * Np_wie * N3;
                        size_t tri_shift_out = i_tri * N_vox;
                                                    
                        for (size_t i_line = 0; i_line < N2; ++i_line) {
                            
                            // create a pointer for quadrature points on current line
                            double * r_src = nullptr;
                            double * t   = nullptr;
                            double L_segment = 0.0;
                            
                            // if current line semgnet is the first line segment in triangle basis
                            if(0 == i_line) {
                                r_src = r_src_p;
                                t     = f_wie_p;
                                L_segment = L_p[i_tri];
                            }
                            else {
                                r_src = r_src_n;
                                t     = f_wie_n;
                                L_segment = L_n[i_tri];
                            }
                            
                            // loop over wire quadrature points
                            for (size_t i_wie = 0; i_wie < Np_wie; ++i_wie) {
                                
                                // compute the radius-vector between the source and observer
                                R_vec[0] = r_obs_voxel[0] - r_src[0 + i_wie * N3 + tri_shift_in];
                                R_vec[1] = r_obs_voxel[1] - r_src[1 + i_wie * N3 + tri_shift_in];
                                R_vec[2] = r_obs_voxel[2] - r_src[2 + i_wie * N3 + tri_shift_in];
                                
                                //compute the distance
                                const double R_dist = vec_norm_l2(R_vec);
                                
                                //get current weights
                                const double weight = vie_quads[i_vie + 1] * vie_quads[j_vie + 1] * vie_quads[k_vie + 1] * wie_quads[i_wie + 1] / 8.0;
                                
                                const double R2   = R_dist * R_dist;
                                const double R3   = R_dist * R2;
                                const double koR  = ko * R_dist;
                                const double koR2 = koR * koR;
                                
                                const dcomplex kappa = exp(-dcomplex(0.0,1.0) * koR) / R3;
                                const dcomplex gamma = dcomplex(0.0, 1.0) * koR + 1.0;
                                const dcomplex P = (koR2 - 3.0 * dcomplex(0.0, 1.0) * koR - 3.0) / R2;
                                const dcomplex Q = (koR2 -       dcomplex(0.0, 1.0) * koR - 1.0) / P;
                                const dcomplex mul_ct = kappa * P;
                                
                                // define dyadic Green's function for K operator
                                const dcomplex  Gxx = R_vec[0] * R_vec[0] - Q;
                                const dcomplex  Gxy = R_vec[0] * R_vec[1];
                                const dcomplex  Gxz = R_vec[0] * R_vec[2];
                                //const dcomplex  Gyy = R_vec[1] * R_vec[1] - Q;
                                //const dcomplex  Gyz = R_vec[1] * R_vec[2];
                                //const dcomplex  Gzz = R_vec[2] * R_vec[2] - Q;
                                //const dcomplex G = kappa * gamma * m_scaling;
                                
                                dcomplex  Kernel_Nop;
                                
                                Kernel_Nop = Gxx * t[0 + i_wie * N3 + tri_shift_in] + Gxy * t[1 + i_wie * N3 + tri_shift_in] + Gxz * t[2 + i_wie * N3 + tri_shift_in]; // Nx
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 1.0 * Kernel_Nop * L_segment/2.0 * em_scaling; // Nx
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[i_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Nx1
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[j_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Nx2
                                dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[k_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Nx3
                                
                                //Kernel_Nop = Gxy * t[0 + i_wie * N3 + tri_shift_in] + Gyy * t[1 + i_wie * N3 + tri_shift_in] + Gyz * t[2 + i_wie * N3 + tri_shift_in]; // Ny
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 1.0 * Kernel_Nop * L_segment/2.0 * em_scaling; // Ny
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[i_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Ny1
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[j_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Ny2
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[k_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Ny3
                                
                                //Kernel_Nop = Gxz * t[0 + i_wie * N3 + tri_shift_in] + Gyz * t[1 + i_wie * N3 + tri_shift_in] + Gzz * t[2 + i_wie * N3 + tri_shift_in]; // Nz
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 1.0 * Kernel_Nop * L_segment/2.0 * em_scaling; // Nz
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[i_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Nz1
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[j_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Nz2
                                //dcomplex update_Zbc_Nop = weight * mul_ct * 0.5*vie_quads[k_vie + Np_vie + 1] * Kernel_Nop * L_segment/2.0 * em_scaling; // Nz3
                                
                                //Kernel_Nop = G * (R_vec[1] * t[2 + i_wie * N3 + tri_shift_in] - R_vec[2] * t[1 + i_wie * N3 + tri_shift_in]); // Kx
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 1.0 * Kernel_Nop; // Kx
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[i_vie + Np_vie + 1] * Kernel_Nop; // Kx1
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[j_vie + Np_vie + 1] * Kernel_Nop; // Kx2
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[k_vie + Np_vie + 1] * Kernel_Nop; // Kx3
                                
                                //Kernel_Nop = G * (R_vec[2] * t[0 + i_wie * N3 + tri_shift_in] - R_vec[0] * t[2 + i_wie * N3 + tri_shift_in]); // Ky
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 1.0 * Kernel_Nop; // Ky
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[i_vie + Np_vie + 1] * Kernel_Nop; // Ky1
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[j_vie + Np_vie + 1] * Kernel_Nop; // Ky2
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[k_vie + Np_vie + 1] * Kernel_Nop; // Ky3
                                
                                //Kernel_Nop = G * (R_vec[0] * t[1 + i_wie * N3 + tri_shift_in] - R_vec[1] * t[0 + i_wie * N3 + tri_shift_in]); // Kz
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 1.0 * Kernel_Nop; // Kz
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[i_vie + Np_vie + 1] * Kernel_Nop; // Kz1
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[j_vie + Np_vie + 1] * Kernel_Nop; // Kz2
                                //dcomplex update_Zbc_Nop = weight * L_segment/2.0 * 0.5*vie_quads[k_vie + Np_vie + 1] * Kernel_Nop; // Kz3
                                
                                // store components
                                Zbc_Nop[i_vox + tri_shift_out].real  += real(update_Zbc_Nop);
                                Zbc_Nop[i_vox + tri_shift_out].imag  += imag(update_Zbc_Nop);
                                
                                
                            }
                        }
                    
                    }
                    
                }
            }
        }
        
    }
    
    // de-allocate memory
    delete [] r_src_p;
    delete [] r_src_n;
    delete [] f_wie_p;
    delete [] f_wie_n;
    delete [] L_p;
    delete [] L_n;
    
}

///////////////////////////////////////////////////////////////////////////
void get_source_coords_mat(double * r_v, double * r_2, double * wie_quads, double * z_wie, double * t, const size_t Np_wie, const size_t N_tri) {
    
    const size_t N3 = 3;
    double lp;
    
    for (size_t i_tri = 0; i_tri < N_tri; ++i_tri) {
        
        // Compute segment vector and length
        double dx = r_v[0 + i_tri*N3] - r_2[0 + i_tri*N3];
        double dy = r_v[1 + i_tri*N3] - r_2[1 + i_tri*N3];
        double dz = r_v[2 + i_tri*N3] - r_2[2 + i_tri*N3];
        double L = sqrt(dx*dx + dy*dy + dz*dz);
        
        // Normalize and apply sign to tangent
        double tx = (dx / L);
        double ty = (dy / L);
        double tz = (dz / L);
        
        for (size_t i_wie = 0; i_wie < Np_wie; ++i_wie) {
            
            // Compute source points
            lp = L/2.0 *(wie_quads[i_wie + Np_wie + 1]+1);

            z_wie[0 + i_wie*N3 + i_tri*Np_wie*N3] = r_2[0 + i_tri*N3] + lp/L * dx;
            z_wie[1 + i_wie*N3 + i_tri*Np_wie*N3] = r_2[1 + i_tri*N3] + lp/L * dy;
            z_wie[2 + i_wie*N3 + i_tri*Np_wie*N3] = r_2[2 + i_tri*N3] + lp/L * dz;
            
            // Set tangent vector (constant for the segment)
            t[0 + i_wie*N3 + i_tri*Np_wie*N3] = lp/L * tx;
            t[1 + i_wie*N3 + i_tri*Np_wie*N3] = lp/L * ty;
            t[2 + i_wie*N3 + i_tri*Np_wie*N3] = lp/L * tz;
        }
    }
}

///////////////////////////////////////////////////////////////////////////
double * compute_segment_lengths(double * r_start, double * r_end, const size_t N_tri) {
    const size_t N3 = 3;
    double * lengths = new double[N_tri];
    
    for (size_t i = 0; i < N_tri; ++i) {
        double dx = r_end[0 + i*N3] - r_start[0 + i*N3];
        double dy = r_end[1 + i*N3] - r_start[1 + i*N3];
        double dz = r_end[2 + i*N3] - r_start[2 + i*N3];
        lengths[i] = sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    return lengths;
}

///////////////////////////////////////////////////////////////////////////
//! function vec_norm_l2 computes the second norm of the 3D vector r
double vec_norm_l2(double * r) { return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);}
