function [MREDM] = pfft_wire_domain(MREDM)

    inp = MREDM.inputs;
    
    %% Get dimensions
    x_3d = squeeze(MREDM.dimensions.r(:,:,:,1));
    y_3d = squeeze(MREDM.dimensions.r(:,:,:,2));
    z_3d = squeeze(MREDM.dimensions.r(:,:,:,3));
    x    = squeeze(MREDM.dimensions.r(:,1,1,1));
    y    = squeeze(MREDM.dimensions.r(1,:,1,2));
    z    = squeeze(MREDM.dimensions.r(1,1,:,3));
    res  = MREDM.dimensions.res;
    idxS = MREDM.dimensions.idxS;
    mask = MREDM.dimensions.mask;
    xc   = MREDM.WIE.coil.F_point(:,1);
    yc   = MREDM.WIE.coil.F_point(:,2);
    zc   = MREDM.WIE.coil.F_point(:,3);
    
    pFFT_kernel                    = 3;
    MREDM.dimensions.pFFT_kernel   = pFFT_kernel;
    MREDM.dimensions.pFFT_3Dkernel = MREDM.dimensions.pFFT_kernel^3;
    
    N_coil2bound = (pFFT_kernel / 2) - 1;
    N_exp_half   = (pFFT_kernel - 1) / 2;
    d_coil2bound = res * N_coil2bound;

    %% Extend the domain
    [r_ext, ext_mask] = pfft_extend_vie_domain(x_3d, y_3d, z_3d, x, y, z, res, idxS, mask, xc, yc, zc, d_coil2bound, N_exp_half);

    %% Store pFFT dimensions
    MREDM.dimensions.pfft_r = r_ext;
    MREDM.dimensions.pfft_n1 = size(r_ext,1);
    MREDM.dimensions.pfft_n2 = size(r_ext,2);
    MREDM.dimensions.pfft_n3 = size(r_ext,3);
    MREDM.dimensions.pfft_mask = ext_mask;
    MREDM.dimensions.pfft_idxS = find(abs(ext_mask(:)));
        
    ext_nS = length(MREDM.dimensions.pfft_idxS);
    ext_nD = MREDM.dimensions.pfft_n1*MREDM.dimensions.pfft_n2*MREDM.dimensions.pfft_n3;
    if inp.PWX == 1
        ext_idx_solution = zeros(12*ext_nS,1);
        for i = 1:12
            ext_idx_solution((i-1)*ext_nS+1:i*ext_nS) = (i-1)*ext_nD + MREDM.dimensions.pfft_idxS;
        end
        ext_idx_domain = (1:12*ext_nD)';
    else
        ext_idx_solution = zeros(3*ext_nS,1);
        for i = 1:3
            ext_idx_solution((i-1)*ext_nS+1:i*ext_nS) = (i-1)*ext_nD + MREDM.dimensions.pfft_idxS;
        end
        ext_idx_domain = (1:3*ext_nD)';
    end
    MREDM.dimensions.pfft_idx_solution = ext_idx_solution;
    MREDM.dimensions.pfft_idx_domain = ext_idx_domain;

    %% Store additional pFFT dimensions
    Dl = vecnorm(MREDM.WIE.coil.S_point(1:end,:)-MREDM.WIE.coil.F_point(1:end,:),2,2);
    line_bounding_box = ceil(Dl./res);

    d_near_vox   = pFFT_kernel + floor(pFFT_kernel / 2);
    d_near_line  = ceil(mean(line_bounding_box));
    d_near       = ceil(1.6 * max(d_near_vox, d_near_line));

    MREDM.dimensions.pfft_kernel_n1 = 2 * d_near + 2 * floor(pFFT_kernel / 2) + mod(pFFT_kernel, 2);
    MREDM.dimensions.pfft_kernel_n2 = 2 * d_near + 2 * floor(pFFT_kernel / 2) + mod(pFFT_kernel, 2);
    MREDM.dimensions.pfft_kernel_n3 = 2 * d_near + 2 * floor(pFFT_kernel / 2) + mod(pFFT_kernel, 2);

    MREDM.dimensions.pfft_near_dist = d_near;

end
