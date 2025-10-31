function [MREDM] = geo_spherical_basis(MREDM)

    n1  = MREDM.dimensions.n1;
    n2  = MREDM.dimensions.n2;
    n3  = MREDM.dimensions.n3;
    res = MREDM.dimensions.res;
    r   = MREDM.dimensions.r;
    dis = inp.Basis_distance;
    thi = inp.Basis_thickness;
    
    Cnt = [r(ceil(n1/2),1,1,1) r(1,ceil(n2/2),1,2) r(1,1,ceil(n3/2),3)];
    dist_rsvd = dis + thi*res;
    x_dom_length     = r(end,1,1,1) - r(1,1,1,1);
    y_dom_length     = r(1,end,1,2) - r(1,1,1,2);
    z_dom_length     = r(1,1,end,3) - r(1,1,1,3);
    max_rad          = sqrt( ((x_dom_length-Cnt(1))/2)^2 + ((y_dom_length-Cnt(2))/2)^2 + ((z_dom_length-Cnt(3))/2)^2);
    x                = [-max_rad+dist_rsvd max_rad+dist_rsvd];
    nx               = ceil((max(x)-min(x))/2/res);
    ny               = nx;
    nz               = nx;
    Cx               = (max(x)+min(x))/2;
    Cy               = Cx;
    Cz               = Cx;
    x                = Cx-nx*res:res:Cx+(nx)*res;
    y                = Cy-ny*res:res:Cy+(ny)*res;
    z                = Cz-nz*res:res:Cz+(nz)*res;
    r_basis          = grid3d(x,y,z);
    r_shell_in       = max_rad + dis;
    r_shell_ou       = max_rad + dist_rsvd;
    sphere_in        = @(r) ( (r(:,:,:,1) - Cnt(1) ).^2 + ( r(:,:,:,2) - Cnt(2) ).^2 + ( r(:,:,:,3) - Cnt(3) ).^2 >= r_shell_in^2);
    sphere_out       = @(r) ( (r(:,:,:,1) - Cnt(1) ).^2 + ( r(:,:,:,2) - Cnt(2) ).^2 + ( r(:,:,:,3) - Cnt(3) ).^2 <= r_shell_ou^2);
    points_in        = sphere_in(r_basis);
    points_out       = sphere_out(r_basis);
    points_surface3D = points_in & points_out;
    basis_idxS       = find(points_surface3D(:));
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
