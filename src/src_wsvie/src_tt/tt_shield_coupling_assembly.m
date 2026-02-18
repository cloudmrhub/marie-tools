function[Zbc_N,Zbc_K] = tt_shield_coupling_assembly(MREDM)

    k0  = MREDM.emc.k0;
    tol = MREDM.inputs.tol_TT;

    Quad_order_vie                           = MREDM.inputs.Quad_order_vie;
    Quad_order_sie                           = MREDM.inputs.Quad_order_sie;
    [wp_vie, z_vie]                          = gauss_1d(Quad_order_vie);
    [Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);
    VIE_quads                                = [Quad_order_vie; wp_vie; z_vie];
    SIE_quads                                = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];

    dims = MREDM.dimensions;
    n1  = dims.n1;
    n2  = dims.n2;
    n3  = dims.n3;
    n4  = dims.N_shield_sie;
    ql  = dims.ql;
    xd  = dims.r(:,:,:,1);
    yd  = dims.r(:,:,:,2);
    zd  = dims.r(:,:,:,3);
    res = dims.res;

    rp = MREDM.SIE.shield.rp;
    rn = MREDM.SIE.shield.rn;
    r2 = MREDM.SIE.shield.r2;
    r3 = MREDM.SIE.shield.r3;
    
    if ql == 12
        funG{1}  = @(J) assemble_coupling_matrices_voxels_N_x (SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{2}  = @(J) assemble_coupling_matrices_voxels_N_x1(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{3}  = @(J) assemble_coupling_matrices_voxels_N_x2(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{4}  = @(J) assemble_coupling_matrices_voxels_N_x3(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{5}  = @(J) assemble_coupling_matrices_voxels_N_y (SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{6}  = @(J) assemble_coupling_matrices_voxels_N_y1(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{7}  = @(J) assemble_coupling_matrices_voxels_N_y2(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{8}  = @(J) assemble_coupling_matrices_voxels_N_y3(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{9}  = @(J) assemble_coupling_matrices_voxels_N_z (SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{10} = @(J) assemble_coupling_matrices_voxels_N_z1(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{11} = @(J) assemble_coupling_matrices_voxels_N_z2(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{12} = @(J) assemble_coupling_matrices_voxels_N_z3(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{13} = @(J) assemble_coupling_matrices_voxels_K_x (SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{14} = @(J) assemble_coupling_matrices_voxels_K_x1(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{15} = @(J) assemble_coupling_matrices_voxels_K_x2(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{16} = @(J) assemble_coupling_matrices_voxels_K_x3(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{17} = @(J) assemble_coupling_matrices_voxels_K_y (SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{18} = @(J) assemble_coupling_matrices_voxels_K_y1(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{19} = @(J) assemble_coupling_matrices_voxels_K_y2(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{20} = @(J) assemble_coupling_matrices_voxels_K_y3(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{21} = @(J) assemble_coupling_matrices_voxels_K_z (SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{22} = @(J) assemble_coupling_matrices_voxels_K_z1(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{23} = @(J) assemble_coupling_matrices_voxels_K_z2(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{24} = @(J) assemble_coupling_matrices_voxels_K_z3(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        print_comp = ["Nx0","Nx1","Nx2","Nx3","Ny0","Ny1","Ny2","Ny3","Nz0","Nz1","Nz2","Nz3",...
                      "Kx0","Kx1","Kx2","Kx3","Ky0","Ky1","Ky2","Ky3","Kz0","Kz1","Kz2","Kz3"];
    else
        funG{1} = @(J) assemble_coupling_matrices_voxels_N_x(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{2} = @(J) assemble_coupling_matrices_voxels_N_y(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{3} = @(J) assemble_coupling_matrices_voxels_N_z(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{4} = @(J) assemble_coupling_matrices_voxels_K_x(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{5} = @(J) assemble_coupling_matrices_voxels_K_y(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        funG{6} = @(J) assemble_coupling_matrices_voxels_K_z(SIE_quads, VIE_quads, k0, n1, n2, n3, xd, yd, zd, res, rp, rn ,r2, r3, J);
        print_comp = ["Nx","Ny","Nz","Kx","Ky","Kz"];
    end
    
    cTT = struct('TT',repmat({[]},2*ql));
    for ij = 1:2*ql
        fprintf('\t\t\tAssembling PWC-%s/Shield',print_comp(ij));
        Green = dmrg_cross_gpu(4,[n1;n2;n3;n4], funG{ij}, tol);
        cTT(ij).TT.d    = Green.d;
        cTT(ij).TT.n    = Green.n;
        cTT(ij).TT.r    = Green.r;
        cTT(ij).TT.ps   = Green.ps;
        cTT(ij).TT.core = Green.core;
        fprintf('\n');
    end
    Zbc_N = cTT(1:ql);
    Zbc_K = cTT(ql+1:2*ql);

end
