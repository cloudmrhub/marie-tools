function [MREDM] = pfft_proj_wire_surface_create_near_lists(MREDM)    

    dims           = MREDM.dimensions;
    wire           = MREDM.WIE.coil;
    coil           = MREDM.SIE.coil;
    N_wie          = dims.N_wie;
    N_sie          = dims.N_sie;
    pFFT_3Dkernel  = dims.pFFT_3Dkernel;
    pFFT_kernel    = dims.pFFT_kernel;
    pfft_near_dist = dims.pfft_near_dist;
    pfft_kernel_n1 = dims.pfft_kernel_n1;
    pfft_kernel_n2 = dims.pfft_kernel_n2;
    pfft_kernel_n3 = dims.pfft_kernel_n3;
    pfft_idxS      = dims.pfft_idxS;
    pfft_n1        = dims.pfft_n1;
    pfft_n2        = dims.pfft_n2;
    pfft_n3        = dims.pfft_n3;
    xd             = squeeze(dims.pfft_r(:,1,1,1));
    yd             = squeeze(dims.pfft_r(1,:,1,2));
    zd             = squeeze(dims.pfft_r(1,1,:,3));
    res            = dims.res; 

    % Initialize
    C2Exp_centers          = zeros(3, N_sie);
    Coil_Func_to_PWX_cells = zeros(pFFT_3Dkernel, N_sie);
    C2Vox_near_list        = zeros(pfft_kernel_n1*pfft_kernel_n2*pfft_kernel_n3, N_sie);
    C2C_near_list          = containers.Map;
    C2B_near_list          = containers.Map;
    
    % Find RWG centers 
    RWG_cntr = pfft_proj_find_RWG_centers(coil);
    % Find Triangle centers 
    Triangle_cntr = wire.S_point;
    % Find all centers
    All_cntr = [Triangle_cntr; RWG_cntr];

    for i_all = 1:N_wie+N_sie

        % Find nearest voxel to RWG function center 
        r_sie = All_cntr(i_all,:).';
        [idx_x_ctr, idx_y_ctr, idx_z_ctr] = pfft_proj_find_nearest_voxel(r_sie, xd, yd, zd, res);
        
        % Find the voxels of the extended domain that RWG functions have to expand to. 
        Coil_Func_to_PWX_cells(:,i_all) = pfft_proj_find_expansion_cell(idx_x_ctr, idx_y_ctr,idx_z_ctr,pFFT_kernel,pfft_n1,pfft_n2,pfft_n3);
   
        % Find the indecies of the extended domain per RWG function that the projection has to be applied on
        C2Vox_near_list(:,i_all)  = pfft_proj_get_near_indecies(idx_x_ctr,idx_y_ctr,idx_z_ctr,pfft_kernel_n1,pfft_n1,pfft_n2,pfft_n3);

        % Store centers of voxels per RWG function
        C2Exp_centers(:,i_all)    = [idx_x_ctr, idx_y_ctr, idx_z_ctr];
        
        % Find the indecies of the expansion voxels that intersect with body voxels 
        C2B_near_vox = intersect(C2Vox_near_list(:,i_all), pfft_idxS);
        if nnz(C2B_near_vox) > 0
            C2B_near_list(num2str(i_all)) = C2B_near_vox;
        end
        
    end

    % Find the closest RWG interactions to apply the precorrection 
    for i_all = 1:N_wie+N_sie
        C2C_near_list(num2str(i_all)) = find(max(abs(C2Exp_centers-C2Exp_centers(:,i_all))) < pfft_near_dist);
    end

    exp_vox_1D                   = (pfft_kernel_n1 - pFFT_kernel) / 2 + (1:pFFT_kernel);
    [i_exp_n1,j_exp_n2,k_exp_n3] = meshgrid(exp_vox_1D, exp_vox_1D, exp_vox_1D);
    index_1D                     = reshape(sub2ind([pfft_kernel_n1 pfft_kernel_n2 pfft_kernel_n3], j_exp_n2, i_exp_n1, k_exp_n3), pFFT_3Dkernel,1);

    MREDM.pfft.index_exp2near.S_1d = index_1D;
    MREDM.pfft.index_exp2near.S_3d = [index_1D; index_1D + pfft_kernel_n1*pfft_kernel_n2*pfft_kernel_n3; index_1D + 2 * pfft_kernel_n1*pfft_kernel_n2*pfft_kernel_n3];
    MREDM.pfft.central_cell_id     = C2Exp_centers;
    MREDM.pfft.cell_id             = Coil_Func_to_PWX_cells;
    MREDM.pfft.near_vox_list       = C2Vox_near_list;
    MREDM.pfft.C2C_near_list       = C2C_near_list;
    MREDM.pfft.C2B_near_list       = C2B_near_list;
    
end