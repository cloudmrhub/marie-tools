function R2 = RHBM_downsample(R, factor, varargin)
% R2 = rhbm_downsample(R, factor, 'Name',Value,...)
% Downsample an RHBM struct isotropically by integer FACTOR (e.g., 2 mm -> 4 mm).
%
% Inputs
%   R           RHBM struct with fields:
%               - r (Nx×Ny×Nz×3), epsilon_r, sigma_e, rho, idxS
%   factor      integer >=1
%
% Name-Value options
%   'UseNearest'  (false)  true = nearest/decimation (keeps sharp labels)
%                           false = block average (physically sensible for eps/sigma/rho)
%   'MaskPolicy'  ('any')  'any' | 'all'  (block reduction rule for idxS)
%   'PreserveFOV' ('crop') 'crop' | 'pad'  (make size divisible by factor)
%   'Verbose'     (false)  print decisions
%
% Output
%   R2          downsampled RHBM (fields: r, epsilon_r, sigma_e, rho, idxS, name)
%
% Notes
% - r is **recomputed** on a coarse rectilinear lattice using dx,dy,dz inferred
%   from R.r. This guarantees spacing = factor × original spacing.
% - 'crop' removes trailing slices so dims are divisible by FACTOR.
%   'pad' replicates the last slice along each axis to keep FOV untouched.
%
% Example
%   R4mm = rhbm_downsample(R2mm, 2, 'PreserveFOV','crop', 'UseNearest',false, 'MaskPolicy','any');

    %% ---- Parse & validate
    p = inputParser;
    p.addRequired('R', @isstruct);
    p.addRequired('factor', @(x) isnumeric(x)&&isscalar(x)&&x>=1&&abs(x-round(x))<eps);
    p.addParameter('UseNearest', false, @(x)islogical(x)&&isscalar(x));
    p.addParameter('MaskPolicy', 'any', @(s) any(strcmpi(s,{'any','all'})));
    p.addParameter('PreserveFOV', 'crop', @(s) any(strcmpi(s,{'crop','pad'})));
    p.addParameter('Verbose', false, @(x)islogical(x)&&isscalar(x));
    p.parse(R, factor, varargin{:});
    useNearest  = p.Results.UseNearest;
    maskPolicy  = lower(p.Results.MaskPolicy);
    preserveFOV = lower(p.Results.PreserveFOV);
    verbose     = p.Results.Verbose;

    assert(isfield(R,'r') && ndims(R.r)==4 && size(R.r,4)==3, 'R.r must be Nx×Ny×Nz×3');
    req = {'epsilon_r','sigma_e','rho','idxS'};
    for k=1:numel(req), assert(isfield(R,req{k}), 'Missing field: %s', req{k}); end

    [Nx,Ny,Nz,~] = size(R.r);

    %% ---- Build mask from idxS and prep arrays
    mask0 = false(Nx,Ny,Nz);
    mask0(R.idxS) = true;

    Aeps  = R.epsilon_r;
    Asig  = R.sigma_e;
    Arho  = R.rho;

    % Decide crop/pad sizes so each dim is divisible by factor
    switch preserveFOV
        case 'crop'
            Nx_e = floor(Nx/factor)*factor;
            Ny_e = floor(Ny/factor)*factor;
            Nz_e = floor(Nz/factor)*factor;

            crop = @(A) A(1:Nx_e, 1:Ny_e, 1:Nz_e);
            Aeps = crop(Aeps); Asig = crop(Asig); Arho = crop(Arho); mask = crop(mask0);

        case 'pad'
            [Aeps, px, py, pz] = pad_to_multiple3d(Aeps, factor);
            [Asig, ~ , ~ , ~ ] = pad_to_multiple3d(Asig, factor);
            [Arho, ~ , ~ , ~ ] = pad_to_multiple3d(Arho, factor);
            [mask, ~ , ~ , ~ ] = pad_to_multiple3d(mask0, factor);
            Nx_e = size(Aeps,1); Ny_e = size(Aeps,2); Nz_e = size(Aeps,3);
            if verbose
                fprintf('Padded by [%d %d %d] to [%d %d %d].\n', px,py,pz, Nx_e,Ny_e,Nz_e);
            end
    end

    Nx2 = Nx_e/factor; Ny2 = Ny_e/factor; Nz2 = Nz_e/factor;
    assert(all(mod([Nx_e Ny_e Nz_e], factor)==0), 'Internal error: sizes not divisible by factor.');

    %% ---- Downsample properties & mask
    if useNearest
        take = @(A) A(1:factor:Nx_e, 1:factor:Ny_e, 1:factor:Nz_e);
        eps2 = take(Aeps);
        sig2 = take(Asig);
        rho2 = take(Arho);
        m2   = take(mask);
    else
        eps2 = block_mean(Aeps, factor, Nx2, Ny2, Nz2);
        sig2 = block_mean(Asig, factor, Nx2, Ny2, Nz2);
        rho2 = block_mean(Arho, factor, Nx2, Ny2, Nz2);
        switch maskPolicy
            case 'any', m2 = block_any(mask, factor, Nx2, Ny2, Nz2);
            case 'all', m2 = block_all(mask, factor, Nx2, Ny2, Nz2);
        end
    end

    %% ---- Recompute r to guarantee exact VIR = factor × original
    % Infer spacings from original rectilinear grid
    x_line = R.r(:,1,1,1);  dx = median(diff(x_line));
    y_line = R.r(1,:,1,2);  dy = median(diff(y_line));
    z_line = R.r(1,1,:,3);  dz = median(diff(z_line));
    x0 = x_line(1); y0 = y_line(1); z0 = z_line(1);

    % New centers
    x2 = x0 + (0:Nx2-1) * (factor*dx);
    y2 = y0 + (0:Ny2-1) * (factor*dy);
    z2 = z0 + (0:Nz2-1) * (factor*dz);
    [X2,Y2,Z2] = ndgrid(x2,y2,z2);

    % Sanity: spacing matches exactly
    tol = 1e-15;
    assert(abs((x2(2)-x2(1)) - factor*dx) < tol ...
        && abs((y2(2)-y2(1)) - factor*dy) < tol ...
        && abs((z2(2)-z2(1)) - factor*dz) < tol, ...
        'Non-uniform spacing detected in R.r (cannot guarantee VIR).');

    %% ---- Assemble output
    R2 = struct();
    R2.r = zeros(Nx2,Ny2,Nz2,3, class(R.r));
    R2.r(:,:,:,1) = cast(X2, class(R.r));
    R2.r(:,:,:,2) = cast(Y2, class(R.r));
    R2.r(:,:,:,3) = cast(Z2, class(R.r));

    R2.epsilon_r = cast(eps2, class(R.epsilon_r));
    R2.sigma_e   = cast(sig2, class(R.sigma_e));
    R2.rho       = cast(rho2, class(R.rho));
    R2.idxS      = find(m2);

    if isfield(R,'name')
        R2.name = sprintf('%s_ds%d', R.name, factor);
    else
        R2.name = sprintf('RHBM_ds%d', factor);
    end
end

%% ----- helpers -----

function B = block_mean(A,f,Nx2,Ny2,Nz2)
    A = reshape(A, f, Nx2, f, Ny2, f, Nz2);
    B = squeeze(mean(mean(mean(A,1,'omitnan'),3,'omitnan'),5,'omitnan'));
end

function B = block_any(A,f,Nx2,Ny2,Nz2)
    A = reshape(A, f, Nx2, f, Ny2, f, Nz2);
    B = squeeze(any(any(any(A,1),3),5));
end

function B = block_all(A,f,Nx2,Ny2,Nz2)
    A = reshape(A, f, Nx2, f, Ny2, f, Nz2);
    B = squeeze(all(all(all(A,1),3),5));
end

function [A2,px,py,pz] = pad_to_multiple3d(A,f)
% Edge-replicate to make each dim a multiple of f (no toolboxes).
    [Nx,Ny,Nz] = size(A);
    rx = mod(f - mod(Nx,f), f);
    ry = mod(f - mod(Ny,f), f);
    rz = mod(f - mod(Nz,f), f);

    A2 = A;
    if rx>0
        A2 = cat(1, A2, repmat(A2(end,:,:), [rx,1,1]));
    end
    if ry>0
        A2 = cat(2, A2, repmat(A2(:,end,:), [1,ry,1]));
    end
    if rz>0
        A2 = cat(3, A2, repmat(A2(:,:,end), [1,1,rz]));
    end
    px=rx; py=ry; pz=rz;
end
