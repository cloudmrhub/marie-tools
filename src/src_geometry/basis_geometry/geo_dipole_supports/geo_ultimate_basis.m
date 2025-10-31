function [MREDM] = geo_ultimate_basis(MREDM)

    inp = MREDM.inputs;
    n1  = MREDM.dimensions.n1;
    n2  = MREDM.dimensions.n2;
    n3  = MREDM.dimensions.n3;
    res = MREDM.dimensions.res;
    dis = inp.Basis_distance;
    thi = inp.Basis_thickness;

    dis = floor(dis/res);
    mask = MREDM.dimensions.mask;
    [size_x, size_y, size_z] = size(mask);
    new_size_x = size_x + 2*dis + 2*thi;
    new_size_y = size_y + 2*dis + 2*thi;
    new_size_z = size_z + 2*dis + 2*thi;
    new_mask = zeros(new_size_x, new_size_y, new_size_z);
    center_x = ceil(new_size_x/2);
    center_y = ceil(new_size_y/2);
    center_z = ceil(new_size_z/2);
    start_x = center_x - floor(size_x/2)+1;
    end_x = center_x + ceil(size_x/2);
    start_y = center_y - floor(size_y/2)+1;
    end_y = center_y + ceil(size_y/2);
    start_z = center_z - floor(size_z/2);
    end_z = center_z + ceil(size_z/2)-1;
    new_mask(start_x:end_x, start_y:end_y, start_z:end_z) = mask;
    disA = imdilate(new_mask, strel('sphere', dis));
    disThiA = imdilate(new_mask, strel('sphere', dis+thi));
    maskB = imsubtract(disThiA, disA);
    maskB = bwareaopen(maskB, 1);
    x = [-new_size_x/2*res new_size_x/2*res];
    y = [-new_size_y/2*res new_size_y/2*res];
    z = [-new_size_z/2*res new_size_z/2*res];
    nx = ceil((max(x)-min(x))/2/res);
    ny = ceil((max(y)-min(y))/2/res);
    nz = ceil((max(z)-min(z))/2/res);
    Cx = (max(x)+min(x))/2;
    Cy = (max(y)+min(y))/2;
    Cz = (max(z)+min(z))/2;
    a = 1; b = 1; c = 1;
    if new_size_x ~= nx*2
        a = 2;
    end
    if new_size_y ~= ny*2
        b = 2;
    end
    if new_size_z ~= nz*2
        c = 2;
    end
    x = Cx-nx*res:res:Cx+(nx-a)*res;
    y = Cy-ny*res:res:Cy+(ny-b)*res;
    z = Cz-nz*res:res:Cz+(nz-c)*res;
    r_basis = grid3d(x,y,z);
    basis_idxS = find(maskB>0);
    
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