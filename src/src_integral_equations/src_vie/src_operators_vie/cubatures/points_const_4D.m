function [X,Y,Z,Xp,Yp,Zp] = points_const_4D(A,B,C,D,r_m,r_n,dx,n,np)

%% inputs and outputs 

% A,B,C,D        : (Np^4 x 1) double vectors. A,B,C,D conatin the points according to Gauss/Clenashaw-Curtis quadrature rule
% r_m            : (3x1) centre of the observation voxel
% r_n            : (3x1) centre of the source voxel
% dx             : resolution
% n              : (3x1) normal vector on the surface of the observation voxel
% np             : (3x1) normal vector on the surface of the source voxel

% X,Y,Z,Xp,Yp,Zp : N^4 integration points for the cases of constant planes {X,Y,Z} U {X',Y',Z'} transformed  for integration over a unit square-square


%% calculate the integration points

% constant observation surface
plane   = find(n ~=0);
% constant source surface
plane_p = find(np~=0);
% length of 4D points = Np^4
N = length(A);

%------------------------------------------

% constant X and constant X'
if plane == 1 && plane_p == 1
    X  = r_m(1) * ones(N,1); 
    Y  = r_m(2) + dx/2*A;
    Z  = r_m(3) + dx/2*B;
    
    Xp = r_n(1) * ones(N,1);
    Yp = r_n(2) + dx/2*C;
    Zp = r_n(3) + dx/2*D;
end

% constant X and constant Y'
if plane == 1 && plane_p == 2
    X  = r_m(1) * ones(N,1); 
    Y  = r_m(2) + dx/2*A;
    Z  = r_m(3) + dx/2*B;
    
    Xp = r_n(1) + dx/2*C;
    Yp = r_n(2) * ones(N,1);
    Zp = r_n(3) + dx/2*D;
end

% constant X and constant Z'
if plane == 1 && plane_p == 3
    X  = r_m(1) * ones(N,1); 
    Y  = r_m(2) + dx/2*A;
    Z  = r_m(3) + dx/2*B;
    
    Xp = r_n(1) + dx/2*C;
    Yp = r_n(2) + dx/2*D;
    Zp = r_n(3) * ones(N,1);
end

%------------------------------------------

% constant Y and constant X'
if plane == 2 && plane_p == 1
    X  = r_m(1) + dx/2*A;
    Y  = r_m(2) * ones(N,1);
    Z  = r_m(3) + dx/2*B;
    
    Xp = r_n(1) * ones(N,1);
    Yp = r_n(2) + dx/2*C;
    Zp = r_n(3) + dx/2*D;
end

% constant Y and constant Y'
if plane == 2 && plane_p == 2
    X  = r_m(1) + dx/2*A;
    Y  = r_m(2) * ones(N,1);
    Z  = r_m(3) + dx/2*B;
    
    Xp = r_n(1) + dx/2*C;
    Yp = r_n(2) * ones(N,1);
    Zp = r_n(3) + dx/2*D;
end

% constant Y and constant Z'
if plane == 2 && plane_p == 3
    X  = r_m(1) + dx/2*A;
    Y  = r_m(2) * ones(N,1);
    Z  = r_m(3) + dx/2*B;
    
    Xp = r_n(1) + dx/2*C;
    Yp = r_n(2) + dx/2*D;
    Zp = r_n(3) * ones(N,1);
end

%------------------------------------------

% constant Z and constant X'
if plane == 3 && plane_p == 1
    X  = r_m(1) + dx/2*A;
    Y  = r_m(2) + dx/2*B;
    Z  = r_m(3) * ones(N,1);
    
    Xp = r_n(1) * ones(N,1);
    Yp = r_n(2) + dx/2*C;
    Zp = r_n(3) + dx/2*D;
end

% constant Z and constant Y'
if plane == 3 && plane_p == 2
    X  = r_m(1) + dx/2*A;
    Y  = r_m(2) + dx/2*B;
    Z  = r_m(3) * ones(N,1);
    
    Xp = r_n(1) + dx/2*C;
    Yp = r_n(2) * ones(N,1);
    Zp = r_n(3) + dx/2*D;
end

% constant Z and constant Z'
if plane == 3 && plane_p == 3
    X  = r_m(1) + dx/2*A;
    Y  = r_m(2) + dx/2*B;
    Z  = r_m(3) * ones(N,1);
    
    Xp = r_n(1) + dx/2*C;
    Yp = r_n(2) + dx/2*D;
    Zp = r_n(3) * ones(N,1);
end

end




