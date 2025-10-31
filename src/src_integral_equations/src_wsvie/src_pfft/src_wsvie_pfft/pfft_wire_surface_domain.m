function [MREDM] = pfft_wire_surface_domain(MREDM)

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
    xc_s = MREDM.SIE.coil.node(1,:);
    yc_s = MREDM.SIE.coil.node(2,:);
    zc_s = MREDM.SIE.coil.node(3,:);
    xc_w = MREDM.WIE.coil.F_point(:,1);
    yc_w = MREDM.WIE.coil.F_point(:,2);
    zc_w = MREDM.WIE.coil.F_point(:,3);
    xc   = [xc_s(:);xc_w(:)];
    yc   = [yc_s(:);yc_w(:)];
    zc   = [zc_s(:);zc_w(:)];

    % Ln                    = MREDM.SIE.coil.Ln;
    % Ln_all                = Ln(:);
    % sorted_edges          = sort(Ln_all, 'descend');
    % max_combined_length_s = sorted_edges(1) + sorted_edges(2);
    % Dl                    = vecnorm(MREDM.WIE.coil.S_point(1:end,:)-MREDM.WIE.coil.F_point(1:end,:),2,2);
    % sorted_edges          = sort(Dl, 'descend');
    % max_combined_length_w = sorted_edges(1) + sorted_edges(2);
    % voxels_required       = ceil(max([max_combined_length_s max_combined_length_w]) / res);
    % % Force voxel count to be ODD
    % if mod(voxels_required,2) == 0
    %     voxels_required = voxels_required + 1;
    % end
    pFFT_kernel                    = 3;%voxels_required + 2;  
    MREDM.dimensions.pFFT_kernel   = pFFT_kernel;
    MREDM.dimensions.pFFT_3Dkernel = MREDM.dimensions.pFFT_kernel^3;
    
    N_coil2bound = (pFFT_kernel / 2) - 1;
    N_exp_half   = (pFFT_kernel - 1) / 2;
    d_coil2bound = res * N_coil2bound;


    %% Extend the domain
    Xb_bb   = zeros(3,2);
    Xc_bb   = zeros(3,2);
    Xext_bb = zeros(3,2);

    Xb_bb(1,1) = min(x_3d(idxS))-res; Xb_bb(1,2) = max(x_3d(idxS))+res;
    Xb_bb(2,1) = min(y_3d(idxS))-res; Xb_bb(2,2) = max(y_3d(idxS))+res;
    Xb_bb(3,1) = min(z_3d(idxS))-res; Xb_bb(3,2) = max(z_3d(idxS))+res;

    Xc_bb(1,1) = min(xc);  
    Xc_bb(1,2) = max(xc);
    Xc_bb(2,1) = min(yc);  
    Xc_bb(2,2) = max(yc);
    Xc_bb(3,1) = min(zc);  
    Xc_bb(3,2) = max(zc);

    %% Define a bounding box for the SCOIL + RHBM region
    for i = 1:size(Xext_bb,1)
        if(Xb_bb(i,1) < Xc_bb(i,1) - d_coil2bound)
            Xext_bb(i,1) = Xb_bb(i,1) - eps;
        else
            n_shift = round(abs(Xb_bb(i,1) - Xc_bb(i,1)) / res);
            Xext_bb(i,1) = Xb_bb(i,1) - res * (n_shift + N_exp_half) - eps;
        end

        if(Xb_bb(i,2) > Xc_bb(i,2) + d_coil2bound)
            Xext_bb(i,2) = Xb_bb(i,2) + eps;
        else
            n_shift = round(abs(Xb_bb(i,2) - Xc_bb(i,2)) / res);
            Xext_bb(i,2) = Xb_bb(i,2) + res * (n_shift + N_exp_half) + eps;
        end
    end
    
    x_ext = Xext_bb(1,1):res:Xext_bb(1,2);
    y_ext = Xext_bb(2,1):res:Xext_bb(2,2);
    z_ext = Xext_bb(3,1):res:Xext_bb(3,2);
    r_ext = grid3d(x_ext,y_ext,z_ext);
    
    %% Find extended domain indexing
    xb_min_new = max(Xb_bb(1,1),min(x)) - res / 3;
    xb_max_new = min(Xb_bb(1,2),max(x)) + res / 3;
    yb_min_new = max(Xb_bb(2,1),min(y)) - res / 3;
    yb_max_new = min(Xb_bb(2,2),max(y)) + res / 3;
    zb_min_new = max(Xb_bb(3,1),min(z)) - res / 3;
    zb_max_new = min(Xb_bb(3,2),max(z)) + res / 3;
    
    idx_vie_x = (x_ext >= xb_min_new) & (x_ext <= xb_max_new);
    idx_vie_y = (y_ext >= yb_min_new) & (y_ext <= yb_max_new);
    idx_vie_z = (z_ext >= zb_min_new) & (z_ext <= zb_max_new);
    
    ext_mask = zeros(size(r_ext,1),size(r_ext,2),size(r_ext,3));
    ext_mask(idx_vie_x,idx_vie_y,idx_vie_z) = mask;

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

    N_sie = MREDM.dimensions.N_sie;
    triag_bounding_box = zeros(N_sie,1);
    for i = 1:N_sie
        phys_line = find(MREDM.SIE.coil.index == i);
        [~, parent_triag] = find(abs(MREDM.SIE.coil.etod) == phys_line);
        triag_bounding_box(i) = ceil(max(MREDM.SIE.coil.Ln(:,parent_triag),[],"all")/res);
    end

    mean_triag_bounding_box = ceil(mean(triag_bounding_box));
    mean_line_bounding_box  = ceil(mean(line_bounding_box));
    d_near_geometry = max(mean_triag_bounding_box, mean_line_bounding_box);

    d_near_vox   = pFFT_kernel + floor(pFFT_kernel / 2);
    d_near       = ceil(1.6 * max(d_near_vox, d_near_geometry));

    MREDM.dimensions.pfft_kernel_n1 = 2 * d_near + 2 * floor(pFFT_kernel / 2) + mod(pFFT_kernel, 2);
    MREDM.dimensions.pfft_kernel_n2 = 2 * d_near + 2 * floor(pFFT_kernel / 2) + mod(pFFT_kernel, 2);
    MREDM.dimensions.pfft_kernel_n3 = 2 * d_near + 2 * floor(pFFT_kernel / 2) + mod(pFFT_kernel, 2);

    MREDM.dimensions.pfft_near_dist = d_near;

end
