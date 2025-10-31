function [BASIS] = sSVD_rwgBasis(MREDM)

    dims           = MREDM.dimensions;
    tol_sSVD       = MREDM.inputs.tol_sSVD;
    tol_DEIM       = MREDM.inputs.tol_DEIM;
    r              = dims.r;
    Nscat          = dims.N_scat;
    idxS           = dims.idxS;
    ql             = dims.ql;
    k0             = MREDM.emc.k0;
    res            = dims.res;
    rp             = MREDM.SIE.coil.rp;
    rn             = MREDM.SIE.coil.rn;
    r2             = MREDM.SIE.coil.r2;
    r3             = MREDM.SIE.coil.r3;
    xds            = r(:,:,:,1);
    yds            = r(:,:,:,2);
    zds            = r(:,:,:,3);
    Quad_order_vie                           = MREDM.inputs.Quad_order_vie_tt;
    Quad_order_sie                           = MREDM.inputs.Quad_order_sie_tt;
    [wp_vie, z_vie]                          = gauss_1d(Quad_order_vie);
    [Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);
    VIE_quads                                = [Quad_order_vie; wp_vie; z_vie];
    SIE_quads                                = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];
    Scoord = [xds(:) yds(:) zds(:)].';

    if ql == 12
        Zbc_Kx  = Assemble_rwg_coupling_matrix_K_x(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kx  = Zbc_Kx(idxS,:);
        Zbc_Kx1 = Assemble_rwg_coupling_matrix_K_x1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kx1 = Zbc_Kx1(idxS,:);
        Zbc_Kx2 = Assemble_rwg_coupling_matrix_K_x2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kx2 = Zbc_Kx2(idxS,:);
        Zbc_Kx3 = Assemble_rwg_coupling_matrix_K_x3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kx3 = Zbc_Kx3(idxS,:);
        Zbc_Ky  = Assemble_rwg_coupling_matrix_K_y(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ky  = Zbc_Ky(idxS,:);
        Zbc_Ky1 = Assemble_rwg_coupling_matrix_K_y1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ky1 = Zbc_Ky1(idxS,:);
        Zbc_Ky2 = Assemble_rwg_coupling_matrix_K_y2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ky2 = Zbc_Ky2(idxS,:);
        Zbc_Ky3 = Assemble_rwg_coupling_matrix_K_y3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ky3 = Zbc_Ky3(idxS,:);
        Zbc_Kz  = Assemble_rwg_coupling_matrix_K_z(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kz  = Zbc_Kz(idxS,:);
        Zbc_Kz1 = Assemble_rwg_coupling_matrix_K_z1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kz1 = Zbc_Kz1(idxS,:);
        Zbc_Kz2 = Assemble_rwg_coupling_matrix_K_z2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kz2 = Zbc_Kz2(idxS,:);
        Zbc_Kz3 = Assemble_rwg_coupling_matrix_K_z3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kz3 = Zbc_Kz3(idxS,:);
        Hinc    = [Zbc_Kx;Zbc_Kx1;Zbc_Kx2;Zbc_Kx3;Zbc_Ky;Zbc_Ky1;Zbc_Ky2;Zbc_Ky3;Zbc_Kz;Zbc_Kz1;Zbc_Kz2;Zbc_Kz3];
        clear Zbc_Kx Zbc_Kx1 Zbc_Kx2 Zbc_Kx3 Zbc_Ky Zbc_Ky1 Zbc_Ky2 Zbc_Ky3 Zbc_Kz Zbc_Kz1 Zbc_Kz2 Zbc_Kz3
        [BASIS.Ub,S_K,V_K] = svd(Hinc,'econ');
        Zbc_Nx  = Assemble_rwg_coupling_matrix_N_x(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nx  = Zbc_Nx(idxS,:);
        Zbc_Nx1 = Assemble_rwg_coupling_matrix_N_x1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nx1 = Zbc_Nx1(idxS,:);
        Zbc_Nx2 = Assemble_rwg_coupling_matrix_N_x2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nx2 = Zbc_Nx2(idxS,:);
        Zbc_Nx3 = Assemble_rwg_coupling_matrix_N_x3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nx3 = Zbc_Nx3(idxS,:);
        Zbc_Ny  = Assemble_rwg_coupling_matrix_N_y(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ny  = Zbc_Ny(idxS,:);
        Zbc_Ny1 = Assemble_rwg_coupling_matrix_N_y1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ny1 = Zbc_Ny1(idxS,:);
        Zbc_Ny2 = Assemble_rwg_coupling_matrix_N_y2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ny2 = Zbc_Ny2(idxS,:);
        Zbc_Ny3 = Assemble_rwg_coupling_matrix_N_y3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ny3 = Zbc_Ny3(idxS,:);
        Zbc_Nz  = Assemble_rwg_coupling_matrix_N_z(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nz  = Zbc_Nz(idxS,:);
        Zbc_Nz1 = Assemble_rwg_coupling_matrix_N_z1(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nz1 = Zbc_Nz1(idxS,:);
        Zbc_Nz2 = Assemble_rwg_coupling_matrix_N_z2(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nz2 = Zbc_Nz2(idxS,:);
        Zbc_Nz3 = Assemble_rwg_coupling_matrix_N_z3(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nz3 = Zbc_Nz3(idxS,:);
        BASIS.Ue = [Zbc_Nx;Zbc_Nx1;Zbc_Nx2;Zbc_Nx3;Zbc_Ny;Zbc_Ny1;Zbc_Ny2;Zbc_Ny3;Zbc_Nz;Zbc_Nz1;Zbc_Nz2;Zbc_Nz3];
        clear Zbc_Nx Zbc_Nx1 Zbc_Nx2 Zbc_Nx3 Zbc_Ny Zbc_Ny1 Zbc_Ny2 Zbc_Ny3 Zbc_Nz Zbc_Nz1 Zbc_Nz2 Zbc_Nz3
    else
        Zbc_Kx  = Assemble_rwg_coupling_matrix_K_x(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kx  = Zbc_Kx(idxS,:);
        Zbc_Ky  = Assemble_rwg_coupling_matrix_K_y(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ky  = Zbc_Ky(idxS,:);
        Zbc_Kz  = Assemble_rwg_coupling_matrix_K_z(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Kz  = Zbc_Kz(idxS,:);
        Hinc    = [Zbc_Kx;Zbc_Ky;Zbc_Kz];
        clear Zbc_Kx Zbc_Ky Zbc_Kz
        [BASIS.Ub,S_K,V_K] = svd(Hinc,'econ');
        Zbc_Nx  = Assemble_rwg_coupling_matrix_N_x(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nx  = Zbc_Nx(idxS,:);
        Zbc_Ny  = Assemble_rwg_coupling_matrix_N_y(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Ny  = Zbc_Ny(idxS,:);
        Zbc_Nz  = Assemble_rwg_coupling_matrix_N_z(Scoord,rp,rn,r2,r3,SIE_quads,VIE_quads,res,k0)*res^3;
        Zbc_Nz  = Zbc_Nz(idxS,:);
        BASIS.Ue = [Zbc_Nx;Zbc_Ny;Zbc_Nz];
        clear Zbc_Nx Zbc_Ny Zbc_Nz
    end

    [BASIS.Ub,S_K,V_K]    = compress_SVD(BASIS.Ub,S_K,V_K,tol_DEIM);
    [P,ndeim,xds,yds,zds] = rSVD_deim(BASIS.Ub,ql,Nscat,idxS,r);
    [BASIS.Ub,Sh,Vh]      = compress_SVD(BASIS.Ub,S_K,V_K,tol_sSVD);
    BASIS.Ue              = BASIS.Ue*Vh/Sh;
    X                     = (P.'*BASIS.Ue)\speye(ndeim,ndeim);

    BASIS.SK  = Sh;
    BASIS.VK  = Vh;
    BASIS.P   = P;
    BASIS.X   = X;
    BASIS.xds = xds;
    BASIS.yds = yds;
    BASIS.zds = zds;

    fprintf('\trank=%i, last singular value = %4.7f.\n',size(BASIS.Ub,2),BASIS.SK(end,end)/BASIS.SK(1,1)); 

end
