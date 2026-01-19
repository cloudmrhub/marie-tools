function [MREDM] = geo_ultimate_basis(MREDM)

    inp = MREDM.inputs;
    n1  = MREDM.dimensions.n1;
    n2  = MREDM.dimensions.n2;
    n3  = MREDM.dimensions.n3;
    res = MREDM.dimensions.res;

    dis = inp.Basis_distance;   
    thi = inp.Basis_thickness;  

    dis = floor(dis/res);

    pad = dis + thi;

    mask = MREDM.dimensions.mask;

    r = MREDM.dimensions.r; 

    x0 = squeeze(r(:,1,1,1));   
    y0 = squeeze(r(1,:,1,2));  
    z0 = squeeze(r(1,1,:,3));   

    y0 = y0(:).';
    z0 = z0(:).';

    x = (x0(1) - pad*res) : res : (x0(end) + pad*res);
    y = (y0(1) - pad*res) : res : (y0(end) + pad*res);
    z = (z0(1) - pad*res) : res : (z0(end) + pad*res);

    new_size_x = numel(x);
    new_size_y = numel(y);
    new_size_z = numel(z);

    new_mask = zeros(new_size_x, new_size_y, new_size_z, 'like', mask);

    new_mask( (pad+1):(pad+n1), (pad+1):(pad+n2), (pad+1):(pad+n3) ) = mask;

    disA    = imdilate(new_mask, strel('sphere', dis));
    disThiA = imdilate(new_mask, strel('sphere', dis + thi));

    maskB = imsubtract(disThiA, disA);
    maskB = bwareaopen(maskB, 1);

    r_basis = grid3d(x, y, z);

    basis_idxS = find(maskB > 0);

    [nb1,nb2,nb3,~]                              = size(r_basis);
    mask_basis                                   = zeros(nb1,nb2,nb3);
    s1                                           = ceil((nb1-n1)/2);
    s2                                           = ceil((nb2-n2)/2);
    s3                                           = ceil((nb3-n3)/2);
    nbv                                          = nb1*nb2*nb3;
    nI                                           = length(basis_idxS);
    basis_idx_solution                           = zeros(MREDM.dimensions.ql*nI,1);

    mask_basis(s1+1:n1+s1,s2+1:n2+s2,s3+1:n3+s3) = MREDM.dimensions.mask;

    basis_idxS_solve                             = find(abs(mask_basis)>0);
    nbI                                          = length(basis_idxS_solve);
    basis_idx_solution_solve                     = zeros(MREDM.dimensions.ql*nbI,1);

    for ii = 1:MREDM.dimensions.ql
        basis_idx_solution( (ii-1)*nI+1:ii*nI ) = (ii-1)*nbv + basis_idxS;
    end
    for ii = 1:MREDM.dimensions.ql
        basis_idx_solution_solve( (ii-1)*nbI+1:ii*nbI ) = (ii-1)*nbv + basis_idxS_solve;
    end

    MREDM.dimensions.mask_basis               = mask_basis;
    MREDM.dimensions.basis_idxS_solve         = basis_idxS_solve;
    MREDM.dimensions.basis_idx_solution_solve = basis_idx_solution_solve;
    MREDM.dimensions.r_basis                  = r_basis;
    MREDM.dimensions.nb1                      = nb1;
    MREDM.dimensions.nb2                      = nb2;
    MREDM.dimensions.nb3                      = nb3;
    MREDM.dimensions.basis_idxS               = basis_idxS;
    MREDM.dimensions.basis_idx_solution       = basis_idx_solution;

end
