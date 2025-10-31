#include <math.h>
#include <complex>
#include <iostream>
#include "directfn.h"
using namespace std;

complex<double> vector_dot(double x[], complex<double> y[])
{
	return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}
/*
complex<double> Kernel  (double rp[],double rq[],Geometry &geom )
{
	double f_g;
	complex<double> Green;

	double k0 = double(2.0)*M_PI;


	f_g = double(1.0);

	    // linear on voxels
    
	double r[3];
    for (int i = 0; i < 3; i++)
	{
		r[i] = rp[i] - rq[i];
	}
	double R = sqrt(vector_dot(r, r));
	// Green's function
	Green = exp(-Iunit * k0 * R) / (double(4.0)*M_PI*R);
	//Green = exp(-Iunit * k0 * R) / (R);
	
	return f_g*Green;
}
complex<double> Kernel  (double rp[],double rq[],Geometry_triangle &geom )
{
	double f_g;
	complex<double> Green;

	double k0 = double(2.0)*M_PI;


	f_g = double(1.0);

	    // linear on voxels
    
	double r[3];
    for (int i = 0; i < 3; i++)
	{
		r[i] = rp[i] - rq[i];
	}
	double R = sqrt(vector_dot(r, r));
	// Green's function
	Green = exp(-Iunit * k0 * R) / (double(4.0)*M_PI*R);
	//Green = exp(-Iunit * k0 * R) / (R);
	
	return f_g*Green;
}
*/
complex<double> Kernel  ( double rp[], double rq[], Geometry &geom)
{
    
	// declaration of variables

    // position vectors
	double rp_c[3], rq_c[3]; // centers of elements
	
    double Np;  // Nodal shape functions
	double Nq;

	double f_g;
	double np[3], nq[3]; // normales
	complex<double> Green;
    complex<double> K;
    
    
	double delta = geom.delta;
    int kernel_type = geom.kerneltype;
	//cout << kernel_type << endl;
    
    
    // rp - rq
    double r[3];
    for (int i = 0; i < 3; i++)
	{
		r[i] = rp[i] - rq[i];
	}
    // |rp - rq|
	double R = sqrt(vector_dot(r, r));
    
	// free space Green's function
	Green = exp(-Iunit * geom.k0 * R) /(double(4.0) * M_PI * R);
	double Green0 = double(1.0)/(double(4.0) * M_PI * R);
	

	if (kernel_type == 0)
		{
			K = Green;    
		}
	else if (kernel_type == 4)
		{
           K = (Green0 - Green)/(geom.k0*geom.k0);
		}
	else 
	{
		// get centers and normals
		for (int i = 0; i < 3; i++)
			{
				rq_c[i] = geom.rq_c[i];
				rp_c[i] = geom.rp_c[i];

				nq[i] = geom.nq[i];
				np[i] = geom.np[i];
			}
		
		if (kernel_type == 1)
		{
            int lp = geom.lp;
            int lq = geom.lq;
			

			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;
			

			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;


            
			f_g = Np*Nq;
            //f_g = Np[lp]*Nq[lq];
			
            K = f_g*Green;
		}   
		else if (kernel_type == 2)
		{
            int lq = geom.lq;
			
            
			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;
            
			f_g = Nq;
			
			//f_g = Nq[lq];
            complex<double> F[3];
			
            for (int i = 0; i < 3; i++)
            {
                F[i] = -r[i]*(- Iunit*geom.k0*Green/R - Green/(R*R) + Green0/(R*R))/(geom.k0*geom.k0);
				
            }
			
            K = f_g*vector_dot(nq, F);
			//cout << K << endl;
			
		}
		else if (kernel_type == 3)
		{
            int lp = geom.lp;
			
			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;
				
            f_g = Np;
			//f_g = Np[lp];
            complex<double> F[3];
			
            for (int i = 0; i < 3; i++)
            {
                F[i] = -r[i]*(- Iunit*geom.k0*Green/R - Green/(R*R) + Green0/(R*R))/(geom.k0*geom.k0);
            }
			
            K = f_g*vector_dot(nq, F);
			
		}
		else if (kernel_type == 5)
		{
            int lp = geom.lp;
            int lq = geom.lq;
			
			
			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;


			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;
            
			f_g = Np*Nq;
            complex<double> F[3];
			
            for (int i = 0; i < 3; i++)
            {
                F[i] = -r[i]*(- Iunit*geom.k0*Green/R - Green/(R*R) + Green0/(R*R))/(geom.k0*geom.k0);
            }
			
            K = f_g*vector_dot(np, F);
		}

		else if (kernel_type == 6)
		{
			
			complex<double> F[3];

			for (int i = 0; i < 3; i++)
			{
				F[i] = -r[i] * (-Iunit*geom.k0*Green / R - Green / (R*R) + Green0 / (R*R)) / (geom.k0*geom.k0) - r[i]*Green0/double(2.0) ;
			}

			K = vector_dot(np, F);

		}

		else if (kernel_type == 7)
		{
            int lp = geom.lp;
			
			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;
				
            f_g = Np;
			K = f_g*(Green0 - Green)/(geom.k0*geom.k0);
		}
		else if (kernel_type == 8)
		{
            int lq = geom.lq;
			
			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;
				
            f_g = Nq;
			K = f_g*(Green0 - Green)/(geom.k0*geom.k0);
		}

	}
	    
    
	
	return K;
}

complex<double> Kernel  ( double rp[], double rq[], Geometry_triangle &geom)
{
    
	// declaration of variables

    // position vectors
	double rp_c[3], rq_c[3]; // centers of elements
	
    double Np;  // Nodal shape functions
	double Nq;

	double f_g;
	double np[3], nq[3]; // normales
	complex<double> Green;
    complex<double> K;
    
    
	double delta = geom.delta;
    int kernel_type = geom.kerneltype;
	
    
    
    // rp - rq
    double r[3];
    for (int i = 0; i < 3; i++)
	{
		r[i] = rp[i] - rq[i];
	}
    // |rp - rq|
	double R = sqrt(vector_dot(r, r));
    
	// free space Green's function
	Green = exp(-Iunit * geom.k0 * R) /(double(4.0) * M_PI * R);
	double Green0 = double(1.0)/(double(4.0) * M_PI * R);

	if (kernel_type == 0)
		{
			K = Green;    
		}
	else if (kernel_type == 4)
		{
           K = (Green0 - Green)/(geom.k0*geom.k0);
		}
	else 
	{
		// get centers end normales
		for (int i = 0; i < 3; i++)
			{
				rq_c[i] = geom.rq_c[i];
				rp_c[i] = geom.rp_c[i];

				nq[i] = geom.nq[i];
				np[i] = geom.np[i];
			}
		
		if (kernel_type == 1)
		{
            int lp = geom.lp;
            int lq = geom.lq;
			
			
			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;


			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;
            
			f_g = Np*Nq;
            //f_g = Np[lp]*Nq[lq];
			
            K = f_g*Green;
		}   
		else if (kernel_type == 2)
		{
            int lq = geom.lq;
			
            
			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;
            
			f_g = Nq;
			//f_g = Nq[lq];
            complex<double> F[3];
			
            for (int i = 0; i < 3; i++)
            {
                F[i] = -r[i]*(- Iunit*geom.k0*Green/R - Green/(R*R) + Green0/(R*R))/(geom.k0*geom.k0);
            }
			
            K = f_g*vector_dot(nq, F);
			
		}
		else if (kernel_type == 3)
		{
            int lp = geom.lp;
			
			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;
				
            f_g = Np;
			//f_g = Np[lp];
            complex<double> F[3];
			
            for (int i = 0; i < 3; i++)
            {
                F[i] = -r[i]*(- Iunit*geom.k0*Green/R - Green/(R*R) + Green0/(R*R))/(geom.k0*geom.k0);
            }
			
            K = f_g*vector_dot(nq, F);
			
		}
		else if (kernel_type == 5)
		{
            int lp = geom.lp;
            int lq = geom.lq;
			
			
			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;


			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;
            
			f_g = Np*Nq;
            complex<double> F[3];
			
            for (int i = 0; i < 3; i++)
            {
                F[i] = -r[i]*(- Iunit*geom.k0*Green/R - Green/(R*R) + Green0/(R*R))/(geom.k0*geom.k0);
            }
			
            K = f_g*vector_dot(np, F);
		} 

		else if (kernel_type == 6)
		{

			complex<double> F[3];

			for (int i = 0; i < 3; i++)
			{
				F[i] = -r[i] * (-Iunit*geom.k0*Green / R - Green / (R*R) + Green0 / (R*R)) / (geom.k0*geom.k0) - r[i] * Green0 / double(2.0);
			}

			K = vector_dot(np, F);

		}

		else if (kernel_type == 7)
		{
            int lp = geom.lp;
			
			if (lp == 0)
				Np = double(1.0);
			else if (lp == 1)
				Np = (rp[0] - rp_c[0]) / delta + double(0.5)*np[0] ;
			else if (lp == 2)
				Np = (rp[1] - rp_c[1]) / delta + double(0.5)*np[1] ;
			else if (lp == 3)
				Np = (rp[2] - rp_c[2]) / delta + double(0.5)*np[2] ;
				
            f_g = Np;
			K = f_g*(Green0 - Green)/(geom.k0*geom.k0);
		}
		else if (kernel_type == 8)
		{
            int lq = geom.lq;
			
			if (lq == 0)
				Nq = double(1.0);
			else if (lq == 1)
				Nq = (rq[0] - rq_c[0]) / delta + double(0.5)*nq[0] ;
			else if (lq == 2)
				Nq = (rq[1] - rq_c[1]) / delta + double(0.5)*nq[1] ;
			else if (lq == 3)
				Nq = (rq[2] - rq_c[2]) / delta + double(0.5)*nq[2] ;
				
            f_g = Nq;
			K = f_g*(Green0 - Green)/(geom.k0*geom.k0);
		}

	}
	    
    
	
	return K;
}