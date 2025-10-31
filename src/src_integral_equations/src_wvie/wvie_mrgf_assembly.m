function[Zbw_N,Zbw_K] = wvie_mrgf_assembly(MREDM)
    
    k0 = MREDM.emc.k0;

    Quad_order_vie = MREDM.inputs.Quad_order_vie;
    Quad_order_wie = MREDM.inputs.Quad_order_wie_coup;
    [wp_vie,z_vie] = gauss_1d(Quad_order_vie);
    [wp_wie,z_wie] = gauss_1d(Quad_order_wie);
    VIE_quads      = [Quad_order_vie; wp_vie; z_vie];
    WIE_quads      = [Quad_order_wie; wp_wie; z_wie];

    xds = h5read(MREDM.inputs.basis_file,'/BASIS/xds');
    yds = h5read(MREDM.inputs.basis_file,'/BASIS/yds');
    zds = h5read(MREDM.inputs.basis_file,'/BASIS/zds');

    dims     = MREDM.dimensions;
    ql       = dims.ql;
    res      = dims.res;
    w_length = length(MREDM.WIE.coil.T_point);
    rp       = MREDM.WIE.coil.F_point(1:w_length,:).';
    r2       = MREDM.WIE.coil.S_point(1:w_length,:).';
    rn       = MREDM.WIE.coil.T_point(1:w_length,:).';
    Scoord   = [xds(:) yds(:) zds(:)].';
    
    if ql == 12
        Zbw_Nx  = Assemble_tri_coupling_matrix_N_x (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nx1 = Assemble_tri_coupling_matrix_N_x1(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nx2 = Assemble_tri_coupling_matrix_N_x2(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nx3 = Assemble_tri_coupling_matrix_N_x3(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ny  = Assemble_tri_coupling_matrix_N_y (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ny1 = Assemble_tri_coupling_matrix_N_y1(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ny2 = Assemble_tri_coupling_matrix_N_y2(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ny3 = Assemble_tri_coupling_matrix_N_y3(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nz  = Assemble_tri_coupling_matrix_N_z (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nz1 = Assemble_tri_coupling_matrix_N_z1(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nz2 = Assemble_tri_coupling_matrix_N_z2(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nz3 = Assemble_tri_coupling_matrix_N_z3(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;    
        Zbw_N   = [Zbw_Nx;Zbw_Nx1;Zbw_Nx2;Zbw_Nx3;Zbw_Ny;Zbw_Ny1;Zbw_Ny2;Zbw_Ny3;Zbw_Nz;Zbw_Nz1;Zbw_Nz2;Zbw_Nz3];
        clear Zbw_Nx Zbw_Nx1 Zbw_Nx2 Zbw_Nx3 Zbw_Ny Zbw_Ny1 Zbw_Ny2 Zbw_Ny3 Zbw_Nz Zbw_Nz1 Zbw_Nz2 Zbw_Nz3
        Zbw_Kx  = Assemble_tri_coupling_matrix_K_x (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kx1 = Assemble_tri_coupling_matrix_K_x1(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kx2 = Assemble_tri_coupling_matrix_K_x2(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kx3 = Assemble_tri_coupling_matrix_K_x3(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ky  = Assemble_tri_coupling_matrix_K_y (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ky1 = Assemble_tri_coupling_matrix_K_y1(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ky2 = Assemble_tri_coupling_matrix_K_y2(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ky3 = Assemble_tri_coupling_matrix_K_y3(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kz  = Assemble_tri_coupling_matrix_K_z (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kz1 = Assemble_tri_coupling_matrix_K_z1(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kz2 = Assemble_tri_coupling_matrix_K_z2(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kz3 = Assemble_tri_coupling_matrix_K_z3(Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_K   = [Zbw_Kx;Zbw_Kx1;Zbw_Kx2;Zbw_Kx3;Zbw_Ky;Zbw_Ky1;Zbw_Ky2;Zbw_Ky3;Zbw_Kz;Zbw_Kz1;Zbw_Kz2;Zbw_Kz3];
        clear Zbw_Kx Zbw_Kx1 Zbw_Kx2 Zbw_Kx3 Zbw_Ky Zbw_Ky1 Zbw_Ky2 Zbw_Ky3 Zbw_Kz Zbw_Kz1 Zbw_Kz2 Zbw_Kz3
    elseif ql == 3
        Zbw_Nx  = Assemble_tri_coupling_matrix_N_x (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ny  = Assemble_tri_coupling_matrix_N_y (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Nz  = Assemble_tri_coupling_matrix_N_z (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_N   = [Zbw_Nx;Zbw_Ny;Zbw_Nz];
        clear Zbw_Nx Zbw_Ny Zbw_Nz
        Zbw_Kx  = Assemble_tri_coupling_matrix_K_x (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Ky  = Assemble_tri_coupling_matrix_K_y (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_Kz  = Assemble_tri_coupling_matrix_K_z (Scoord,rp,rn,r2,WIE_quads,VIE_quads,res,k0) * res^3;
        Zbw_K   = [Zbw_Kx;Zbw_Ky;Zbw_Kz];
        clear Zbw_Kx Zbw_Ky Zbw_Kz
    end

end
