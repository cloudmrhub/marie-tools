function I_VA = singular_VA_lin(l,lp,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4,r5,r6,r7)

%% inputs and outputs

% l                    : \in {1,2,3,4} indicating the scalar part of the linear/pulse testing function for each component of the vector field
% lp                   : \in {1,2,3,4} indicating the scalar part of the linear/pulse basis function for each component of the vector field
% k0                   : wavenumber
% dx                   : resolution
% Np                   : Number of 1D integration points. Also used for the singular integrals
% ker_type             : \in {1,2,3,4} indicating the reduced surface-surface kernel to be calculated
% np                   : (3x1) normal vector on the surface of the observation voxel
% nq                   : (3x1) normal vector on the surface of the source voxel
% r1,r2,r3,r4,r5,r6,r7 : (3x1) ordered points according to DIRECTFN's logic

% I_VA                 : complex double vertex adjacent singular integral


%% example 

% k0 = 0.2*pi;
% Np = 5;
% dx = 1;

% r1 = [0., 0., 0.];
% r2 = [0.,dx,0.];
% r3 = [0., dx, dx];
% r4 = [0., 0., dx];
% r5 = [0., 0., 2*dx];
% r6 = [0., dx, 2*dx];
% nq = [0.,0.,1.];
% np = [0.,0.,1.];
% rq_c = [dx/2,dx/2,0.];
% rp_c = [dx/2, 3*dx/2, 0.];


%% call C++ mex file to calculate Vertex Adjacent singular integral

% centre of squares
rq_c = (r1 + r3)/2;
rp_c = (r3 + r6)/2;

I_VA = solve_va(r1,r2,r3,r4,r5,r6,r7,Np,Np,Np,Np,k0,dx,rq_c,rp_c,nq,np,ker_type,l,lp);

end