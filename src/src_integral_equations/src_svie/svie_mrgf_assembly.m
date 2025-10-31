function[Zbc_N,Zbc_K] = svie_mrgf_assembly(MREDM,coil)
    
    k0 = MREDM.emc.k0;

    Quad_order_vie                           = MREDM.inputs.Quad_order_vie;
    Quad_order_sie                           = MREDM.inputs.Quad_order_sie;
    [wp_vie, z_vie]                          = gauss_1d(Quad_order_vie);
    [Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);
    VIE_quads                                = [Quad_order_vie; wp_vie; z_vie];
    SIE_quads                                = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];

    xds = h5read(MREDM.inputs.basis_file,'/BASIS/xds');
    yds = h5read(MREDM.inputs.basis_file,'/BASIS/yds');
    zds = h5read(MREDM.inputs.basis_file,'/BASIS/zds');
    
    dims   = MREDM.dimensions;
    ql     = dims.ql;
    res    = dims.res;
    rp     = coil.rp;
    rn     = coil.rn;
    r2     = coil.r2;
    r3     = coil.r3;
    Scoord = [xds(:) yds(:) zds(:)].';

    if ql == 12
        Zbc_Nx  = Assemble_rwg_coupling_matrix_N_x (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nx1 = Assemble_rwg_coupling_matrix_N_x1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nx2 = Assemble_rwg_coupling_matrix_N_x2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nx3 = Assemble_rwg_coupling_matrix_N_x3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ny  = Assemble_rwg_coupling_matrix_N_y (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ny1 = Assemble_rwg_coupling_matrix_N_y1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ny2 = Assemble_rwg_coupling_matrix_N_y2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ny3 = Assemble_rwg_coupling_matrix_N_y3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nz  = Assemble_rwg_coupling_matrix_N_z (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nz1 = Assemble_rwg_coupling_matrix_N_z1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nz2 = Assemble_rwg_coupling_matrix_N_z2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nz3 = Assemble_rwg_coupling_matrix_N_z3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_N   = [Zbc_Nx;Zbc_Nx1;Zbc_Nx2;Zbc_Nx3;Zbc_Ny;Zbc_Ny1;Zbc_Ny2;Zbc_Ny3;Zbc_Nz;Zbc_Nz1;Zbc_Nz2;Zbc_Nz3];
        clear Zbc_Nx Zbc_Nx1 Zbc_Nx2 Zbc_Nx3 Zbc_Ny Zbc_Ny1 Zbc_Ny2 Zbc_Ny3 Zbc_Nz Zbc_Nz1 Zbc_Nz2 Zbc_Nz3
        Zbc_Kx  = Assemble_rwg_coupling_matrix_K_x (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kx1 = Assemble_rwg_coupling_matrix_K_x1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kx2 = Assemble_rwg_coupling_matrix_K_x2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kx3 = Assemble_rwg_coupling_matrix_K_x3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ky  = Assemble_rwg_coupling_matrix_K_y (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ky1 = Assemble_rwg_coupling_matrix_K_y1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ky2 = Assemble_rwg_coupling_matrix_K_y2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ky3 = Assemble_rwg_coupling_matrix_K_y3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kz  = Assemble_rwg_coupling_matrix_K_z (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kz1 = Assemble_rwg_coupling_matrix_K_z1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kz2 = Assemble_rwg_coupling_matrix_K_z2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kz3 = Assemble_rwg_coupling_matrix_K_z3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_K   = [Zbc_Kx;Zbc_Kx1;Zbc_Kx2;Zbc_Kx3;Zbc_Ky;Zbc_Ky1;Zbc_Ky2;Zbc_Ky3;Zbc_Kz;Zbc_Kz1;Zbc_Kz2;Zbc_Kz3];
        clear Zbc_Kx Zbc_Kx1 Zbc_Kx2 Zbc_Kx3 Zbc_Ky Zbc_Ky1 Zbc_Ky2 Zbc_Ky3 Zbc_Kz Zbc_Kz1 Zbc_Kz2 Zbc_Kz3
    else
        Zbc_Nx  = Assemble_rwg_coupling_matrix_N_x (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ny  = Assemble_rwg_coupling_matrix_N_y (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Nz  = Assemble_rwg_coupling_matrix_N_z (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_N   = [Zbc_Nx;Zbc_Ny;Zbc_Nz];
        clear Zbc_Nx Zbc_Ny Zbc_Nz
        Zbc_Kx  = Assemble_rwg_coupling_matrix_K_x (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Ky  = Assemble_rwg_coupling_matrix_K_y (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_Kz  = Assemble_rwg_coupling_matrix_K_z (Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0) * res^3;
        Zbc_K   = [Zbc_Kx;Zbc_Ky;Zbc_Kz];
        clear Zbc_Kx Zbc_Ky Zbc_Kz
    end

end
