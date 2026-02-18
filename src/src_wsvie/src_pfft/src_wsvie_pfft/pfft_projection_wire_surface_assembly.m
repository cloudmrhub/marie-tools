function [PS,MREDM] = pfft_projection_wire_surface_assembly(MREDM)
    
    MREDM = pfft_proj_wire_surface_create_near_lists(MREDM);
    coil  = MREDM.SIE.coil;
    wire  = MREDM.WIE.coil;
    k0    = MREDM.emc.k0;
    ce    = MREDM.emc.ce;
    
    dims            = MREDM.dimensions;
    pfft_r          = dims.pfft_r;
    res             = dims.res;
    pfft_idxS       = dims.pfft_idxS;
    col_distance    = (dims.pfft_kernel_n1 - 1) / 2 + 1;
    central_cell_id = MREDM.pfft.central_cell_id;
    cell_id         = MREDM.pfft.cell_id;
    pFFT_kernel     = dims.pFFT_kernel;
    pFFT_3Dkernel   = dims.pFFT_3Dkernel;
    N_sie           = dims.N_sie;
    N_wie           = dims.N_wie;
    Npfft           = dims.pfft_n1*dims.pfft_n2*dims.pfft_n3;
    N_cells         = dims.pFFT_3Dkernel;
    PWX_sz          = dims.ql;
    N_basis         = dims.l;
    Next_dofs       = Npfft * PWX_sz;
    N_vie_dofs      = size(pfft_idxS,1) * PWX_sz;

    % Get points on Lebedev sphere 
    N_colloc_pts   = 26;
    Lebedev_Sphere = getLebedevSphere(N_colloc_pts);
    col_points     = col_distance .* res .* [Lebedev_Sphere.x.'; Lebedev_Sphere.y.'; Lebedev_Sphere.z.'];

    % Get Quadrature Points 
    Quad_order_C2Col                         = 1;
    Quad_order_sie                           = MREDM.inputs.Quad_order_sie;
    Quad_order_wie                           = MREDM.inputs.Quad_order_wie_coup;
    Quad_order_vie                           = MREDM.inputs.Quad_order_vie;
    [wp_vie, z_vie]                          = gauss_1d(Quad_order_C2Col);
    [Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);
    [wp_wie, z_wie]                          = gauss_1d(Quad_order_wie);
    VIE_quads                                = [Quad_order_C2Col; wp_vie; z_vie];
    SIE_quads                                = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];
    WIE_quads                                = [Quad_order_wie; wp_wie; z_wie];

    % Get centers of voxels for the projection
    r_center    = grid3d(res*(linspace(0,pFFT_kernel-1,pFFT_kernel)-floor(pFFT_kernel / 2)));
    vox_centers = round(reshape(r_center,[pFFT_3Dkernel, 3]),4);
  
    % Return the 3x3 Dyadic Green function kernel that models the interactions between collocation points in the Lebedev Sphere and voxels
    A     = pfft_proj_pwx_to_collocation(vox_centers.', res, col_points, Quad_order_vie, PWX_sz, N_basis, k0, ce);
    Val   = cell(N_wie+N_sie,1);
    I_val = cell(N_wie+N_sie,1);

    for i_wiesie = 1:N_wie+N_sie

        ctr_vox_idx     = central_cell_id(:, i_wiesie);
        cel_vox_idx     = cell_id(:, i_wiesie);
        r_c             = squeeze(pfft_r(ctr_vox_idx(1), ctr_vox_idx(2), ctr_vox_idx(3),:));
        cur_col_points  = col_points + r_c; 
        
        if i_wiesie <= N_wie
            % Finds all three nodes associated with a triangle basis function  
            rp = wire.F_point(i_wiesie,:).';
            r2 = wire.S_point(i_wiesie,:).';
            rn = wire.T_point(i_wiesie,:).';
            % Find the rhs for the projection 
            rhs1 = Assemble_tri_coupling_matrix_N_x(cur_col_points,rp,rn,r2,WIE_quads,VIE_quads,res,k0);
            rhs2 = Assemble_tri_coupling_matrix_N_y(cur_col_points,rp,rn,r2,WIE_quads,VIE_quads,res,k0);
            rhs3 = Assemble_tri_coupling_matrix_N_z(cur_col_points,rp,rn,r2,WIE_quads,VIE_quads,res,k0);
            rhs = [rhs1(:);rhs2(:);rhs3(:)];
        else
            % Finds all four vertices associated with an RWG basis function  
            rp = coil.rp(:,i_wiesie-N_wie);
            rn = coil.rn(:,i_wiesie-N_wie);
            r2 = coil.r2(:,i_wiesie-N_wie);
            r3 = coil.r3(:,i_wiesie-N_wie);
            % Find the rhs for the projection 
            rhs1 = Assemble_rwg_coupling_matrix_N_x(cur_col_points,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0);
            rhs2 = Assemble_rwg_coupling_matrix_N_y(cur_col_points,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0);
            rhs3 = Assemble_rwg_coupling_matrix_N_z(cur_col_points,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0);
            rhs = [rhs1(:);rhs2(:);rhs3(:)];
        end
        
        % Solve the system to find the value of currents in the projection voxels
        Ic = A \ rhs;
        
        count_vie = 1;
        j_vie = zeros(N_cells * PWX_sz,1);
        k_vie = zeros(N_cells * PWX_sz,1);

        for q = 1:3
            for basis = 1:N_basis
                for vox_num = 1:N_cells
                    vox_idx = cel_vox_idx(vox_num);
                    j_vie(count_vie) = vox_idx + Npfft * (basis - 1) + Npfft * N_basis * (q - 1);
                    k_vie(count_vie) = vox_num + N_cells * (basis-1) + N_cells * N_basis * (q - 1);
                    count_vie = count_vie + 1;
                end
            end
        end

        Val{i_wiesie}   = Ic(k_vie);
        I_val{i_wiesie} = j_vie;

    end

    Val   = vertcat(Val{:});
    I_val = vertcat(I_val{:});

    P     = sparse(I_val, (1:N_wie+N_sie) .* ones(N_cells * PWX_sz, 1), Val, Next_dofs, N_wie+N_sie);

    I_Jb_idx = repmat(pfft_idxS.',1, PWX_sz) + kron(0:Npfft:Npfft*PWX_sz-Npfft,ones(size(pfft_idxS.')));
    J_Jb_idx = 1:size(I_Jb_idx,2);
    V_Jb_idx = ones(size(I_Jb_idx));
    
    S = sparse(I_Jb_idx, J_Jb_idx, V_Jb_idx, Next_dofs, N_vie_dofs);

    PS = [P S];

end