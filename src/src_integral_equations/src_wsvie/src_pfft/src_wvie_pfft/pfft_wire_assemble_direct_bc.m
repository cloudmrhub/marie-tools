function [Zc2b_N,Zc2b_K] = pfft_wire_assemble_direct_bc(MREDM)

    k0             = MREDM.emc.k0;
    wire           = MREDM.WIE.coil;
    C2B_near_list  = MREDM.pfft.C2B_near_list;
    dims           = MREDM.dimensions;
    pfft_idxS      = dims.pfft_idxS;
    N_b            = dims.N_b;
    N_wie          = dims.N_wie;
    N_scat         = dims.N_scat;
    res            = dims.res;  
    ql             = dims.ql;
    xd             = squeeze(dims.pfft_r(:,:,:,1));
    yd             = squeeze(dims.pfft_r(:,:,:,2));
    zd             = squeeze(dims.pfft_r(:,:,:,3));
    Quad_order_wie = MREDM.inputs.Quad_order_wie_coup;
    Quad_order_vie = MREDM.inputs.Quad_order_vie;

    [wp_vie, z_vie] = gauss_1d(Quad_order_vie);
    [wp_wie, z_wie] = gauss_1d(Quad_order_wie);
    VIE_quads = [Quad_order_vie; wp_vie; z_vie];
    WIE_quads = [Quad_order_wie; wp_wie; z_wie];

    Zbc_Nop_val = cell(C2B_near_list.Count,1);
    Zbc_Kop_val = cell(C2B_near_list.Count,1);
    I_idx       = cell(C2B_near_list.Count,1);
    J_idx       = cell(C2B_near_list.Count,1);

    for key = C2B_near_list.keys
        i_wie                  = str2double(key{:});
        near_scat_vox_ext      = C2B_near_list(key{:});
        Scoord                 = [xd(near_scat_vox_ext), yd(near_scat_vox_ext), zd(near_scat_vox_ext)];
        [~, index_scat2dofs]   = ismember(near_scat_vox_ext,pfft_idxS);
        idx_B_dofs             = repmat(index_scat2dofs,ql,1) + kron(N_scat*(0:ql-1).', ones(size(index_scat2dofs)));
        rp_i                   = wire.F_point(i_wie,:).';
        r2_i                   = wire.S_point(i_wie,:).';
        rn_i                   = wire.T_point(i_wie,:).';
                               
        if ql == 3
            Zbc_Nx  = Assemble_tri_coupling_matrix_N_x (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ny  = Assemble_tri_coupling_matrix_N_y (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nz  = Assemble_tri_coupling_matrix_N_z (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nop_i = [Zbc_Nx;Zbc_Ny;Zbc_Nz];
            Zbc_Kx  = Assemble_tri_coupling_matrix_K_x (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ky  = Assemble_tri_coupling_matrix_K_y (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kz  = Assemble_tri_coupling_matrix_K_z (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kop_i = [Zbc_Kx;Zbc_Ky;Zbc_Kz];
        elseif ql == 12
            Zbc_Nx  = Assemble_tri_coupling_matrix_N_x (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nx1 = Assemble_tri_coupling_matrix_N_x1(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nx2 = Assemble_tri_coupling_matrix_N_x2(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nx3 = Assemble_tri_coupling_matrix_N_x3(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ny  = Assemble_tri_coupling_matrix_N_y (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ny1 = Assemble_tri_coupling_matrix_N_y1(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ny2 = Assemble_tri_coupling_matrix_N_y2(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ny3 = Assemble_tri_coupling_matrix_N_y3(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nz  = Assemble_tri_coupling_matrix_N_z (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nz1 = Assemble_tri_coupling_matrix_N_z1(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nz2 = Assemble_tri_coupling_matrix_N_z2(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nz3 = Assemble_tri_coupling_matrix_N_z3(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Nop_i = [Zbc_Nx;Zbc_Nx1;Zbc_Nx2;Zbc_Nx3;Zbc_Ny;Zbc_Ny1;Zbc_Ny2;Zbc_Ny3;Zbc_Nz;Zbc_Nz1;Zbc_Nz2;Zbc_Nz3];
            Zbc_Kx  = Assemble_tri_coupling_matrix_K_x (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kx1 = Assemble_tri_coupling_matrix_K_x1(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kx2 = Assemble_tri_coupling_matrix_K_x2(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kx3 = Assemble_tri_coupling_matrix_K_x3(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ky  = Assemble_tri_coupling_matrix_K_y (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ky1 = Assemble_tri_coupling_matrix_K_y1(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ky2 = Assemble_tri_coupling_matrix_K_y2(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Ky3 = Assemble_tri_coupling_matrix_K_y3(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kz  = Assemble_tri_coupling_matrix_K_z (Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kz1 = Assemble_tri_coupling_matrix_K_z1(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kz2 = Assemble_tri_coupling_matrix_K_z2(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kz3 = Assemble_tri_coupling_matrix_K_z3(Scoord.',rp_i,rn_i,r2_i,WIE_quads,VIE_quads,res,k0)*res^3;
            Zbc_Kop_i = [Zbc_Kx;Zbc_Kx1;Zbc_Kx2;Zbc_Kx3;Zbc_Ky;Zbc_Ky1;Zbc_Ky2;Zbc_Ky3;Zbc_Kz;Zbc_Kz1;Zbc_Kz2;Zbc_Kz3];
        end
        Zbc_Nop_val{i_wie}     = Zbc_Nop_i;
        Zbc_Kop_val{i_wie}     = Zbc_Kop_i;        
        I_idx{i_wie}           = idx_B_dofs;
        J_idx{i_wie}           = i_wie * ones(size(Zbc_Nop_i));  
    end

    Zbc_Nop_val = vertcat(Zbc_Nop_val{:});
    Zbc_Kop_val = vertcat(Zbc_Kop_val{:});
    I_idx       = vertcat(I_idx{:});
    J_idx       = vertcat(J_idx{:});

    % form sparse coupling matrix
    Zc2b_N = sparse(I_idx, J_idx, Zbc_Nop_val, N_b, N_wie);
    Zc2b_K = sparse(I_idx, J_idx, Zbc_Kop_val, N_b, N_wie);

end