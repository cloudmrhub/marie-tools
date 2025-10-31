function r = grid3d(x, y, z)
%%    Generates a 3D grid
% _________________________________________________________________________
%
%       Generates a 3D cartesian grid or coordinates
%
% _________________________________________________________________________
%
%% INPUT
%   x           positions of the x coordinates
%
%
%% OPTIONAL INPUT
%   y           positions of the y coordinates
%   z           positions of the z coordinates
%
%
%% OUTPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if nargin == 1
    y = x;
    z = x;
end

% define the dimensions
L = length(x);
M = length(y);
N = length(z);

% allocate space
r = zeros(L,M,N,3);

% -------------------------------------------------------------------------
% Fill data
% -------------------------------------------------------------------------

% generate uniform domain
[X,Y,Z] = meshgrid(x,y,z);

% permute tensor to inforce row-major order
X = permute(X, [2, 1, 3, 4]);
Y = permute(Y, [2, 1, 3, 4]);
Z = permute(Z, [2, 1, 3, 4]);

% assing tensors of axis components to r
r(:,:,:,1) = X;
r(:,:,:,2) = Y;
r(:,:,:,3) = Z;


