function[Zbc_N,Zbc_K,rank_min_N,rank_min_K] = svie_coupling_assembly_multicolumn_Tucker(MREDM,coil,N_sie)
    
    k0 = MREDM.emc.k0;

    tol                                      = MREDM.inputs.tol_HOSVD_couple;
    Quad_order_vie                           = MREDM.inputs.Quad_order_vie;
    Quad_order_sie                           = MREDM.inputs.Quad_order_sie;
    [wp_vie, z_vie]                          = gauss_1d(Quad_order_vie);
    [Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);
    VIE_quads                                = [Quad_order_vie; wp_vie; z_vie];
    SIE_quads                                = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];

    dims   = MREDM.dimensions;
    ql     = dims.ql;
    xd     = dims.r(:,:,:,1);
    yd     = dims.r(:,:,:,2);
    zd     = dims.r(:,:,:,3);
    n1     = dims.n1;
    n2     = dims.n2;
    n3     = dims.n3;
    res    = dims.res;
    Scoord = [xd(:) yd(:) zd(:)].';
    
    rp = coil.rp;
    rn = coil.rn;
    r2 = coil.r2;
    r3 = coil.r3;

    Zbc_N = struct('G',repmat({[]},N_sie,ql), 'U1',repmat({[]},N_sie,ql), 'U2',repmat({[]},N_sie,ql), 'U3',repmat({[]},N_sie,ql));
    Zbc_K = struct('G',repmat({[]},N_sie,ql), 'U1',repmat({[]},N_sie,ql), 'U2',repmat({[]},N_sie,ql), 'U3',repmat({[]},N_sie,ql));

    if ql == 12
        for i = 1:N_sie
            [Zbc_N(i,1).G,Zbc_N(i,1).U1,Zbc_N(i,1).U2,Zbc_N(i,1).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_x (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,2).G,Zbc_N(i,2).U1,Zbc_N(i,2).U2,Zbc_N(i,2).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_x1(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,3).G,Zbc_N(i,3).U1,Zbc_N(i,3).U2,Zbc_N(i,3).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_x2(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,4).G,Zbc_N(i,4).U1,Zbc_N(i,4).U2,Zbc_N(i,4).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_x3(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,5).G,Zbc_N(i,5).U1,Zbc_N(i,5).U2,Zbc_N(i,5).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_y (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,6).G,Zbc_N(i,6).U1,Zbc_N(i,6).U2,Zbc_N(i,6).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_y1(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,7).G,Zbc_N(i,7).U1,Zbc_N(i,7).U2,Zbc_N(i,7).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_y2(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,8).G,Zbc_N(i,8).U1,Zbc_N(i,8).U2,Zbc_N(i,8).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_y3(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,9).G,Zbc_N(i,9).U1,Zbc_N(i,9).U2,Zbc_N(i,9).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_z (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,10).G,Zbc_N(i,10).U1,Zbc_N(i,10).U2,Zbc_N(i,10).U3] = hosvd(reshape(Assemble_rwg_coupling_matrix_N_z1(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,11).G,Zbc_N(i,11).U1,Zbc_N(i,11).U2,Zbc_N(i,11).U3] = hosvd(reshape(Assemble_rwg_coupling_matrix_N_z2(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_N(i,12).G,Zbc_N(i,12).U1,Zbc_N(i,12).U2,Zbc_N(i,12).U3] = hosvd(reshape(Assemble_rwg_coupling_matrix_N_z3(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            
            [Zbc_K(i,1).G,Zbc_K(i,1).U1,Zbc_K(i,1).U2,Zbc_K(i,1).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_x (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,2).G,Zbc_K(i,2).U1,Zbc_K(i,2).U2,Zbc_K(i,2).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_x1(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,3).G,Zbc_K(i,3).U1,Zbc_K(i,3).U2,Zbc_K(i,3).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_x2(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,4).G,Zbc_K(i,4).U1,Zbc_K(i,4).U2,Zbc_K(i,4).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_x3(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,5).G,Zbc_K(i,5).U1,Zbc_K(i,5).U2,Zbc_K(i,5).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_y (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,6).G,Zbc_K(i,6).U1,Zbc_K(i,6).U2,Zbc_K(i,6).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_y1(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,7).G,Zbc_K(i,7).U1,Zbc_K(i,7).U2,Zbc_K(i,7).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_y2(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,8).G,Zbc_K(i,8).U1,Zbc_K(i,8).U2,Zbc_K(i,8).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_y3(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,9).G,Zbc_K(i,9).U1,Zbc_K(i,9).U2,Zbc_K(i,9).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_z (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,10).G,Zbc_K(i,10).U1,Zbc_K(i,10).U2,Zbc_K(i,10).U3] = hosvd(reshape(Assemble_rwg_coupling_matrix_K_z1(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,11).G,Zbc_K(i,11).U1,Zbc_K(i,11).U2,Zbc_K(i,11).U3] = hosvd(reshape(Assemble_rwg_coupling_matrix_K_z2(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbc_K(i,12).G,Zbc_K(i,12).U1,Zbc_K(i,12).U2,Zbc_K(i,12).U3] = hosvd(reshape(Assemble_rwg_coupling_matrix_K_z3(Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
        end
    elseif ql == 3
        for i = 1:N_sie
            [Zbc_N(i,1).G,Zbc_N(i,1).U1,Zbc_N(i,1).U2,Zbc_N(i,1).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_x (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbc_N(i,2).G,Zbc_N(i,2).U1,Zbc_N(i,2).U2,Zbc_N(i,2).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_y (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbc_N(i,3).G,Zbc_N(i,3).U1,Zbc_N(i,3).U2,Zbc_N(i,3).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_N_z (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbc_K(i,1).G,Zbc_K(i,1).U1,Zbc_K(i,1).U2,Zbc_K(i,1).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_x (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbc_K(i,2).G,Zbc_K(i,2).U1,Zbc_K(i,2).U2,Zbc_K(i,2).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_y (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbc_K(i,3).G,Zbc_K(i,3).U1,Zbc_K(i,3).U2,Zbc_K(i,3).U3]     = hosvd(reshape(Assemble_rwg_coupling_matrix_K_z (Scoord,rp(:,i),rn(:,i),r2(:,i),r3(:,i),SIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
        end
    end

    % Crop to desired ranks
    rank_min_N = max([n1 n2 n3]);
    rank_min_K = max([n1 n2 n3]);
    for i = 1:N_sie
        for j = 1:ql
            [r1,r2,r3] = size(Zbc_N(i,j).G);
            if min([r1 r2 r3]) < rank_min_N
                rank_min_N = min([r1 r2 r3]);
            end
            [r1,r2,r3] = size(Zbc_K(i,j).G);
            if min([r1 r2 r3]) < rank_min_K
                rank_min_K = min([r1 r2 r3]);
            end
        end
    end
    for i = 1:N_sie
        for j = 1:ql
            Zbc_N(i,j).G = reshape(Zbc_N(i,j).G(1:rank_min_N,1:rank_min_N,1:rank_min_N),rank_min_N,rank_min_N^2);
            Zbc_N(i,j).U1 = Zbc_N(i,j).U1(:,1:rank_min_N);
            Zbc_N(i,j).U2 = Zbc_N(i,j).U2(:,1:rank_min_N);
            Zbc_N(i,j).U3 = Zbc_N(i,j).U3(:,1:rank_min_N).';
            Zbc_K(i,j).G = reshape(Zbc_K(i,j).G(1:rank_min_K,1:rank_min_K,1:rank_min_K),rank_min_K,rank_min_K^2);
            Zbc_K(i,j).U1 = Zbc_K(i,j).U1(:,1:rank_min_K);
            Zbc_K(i,j).U2 = Zbc_K(i,j).U2(:,1:rank_min_K);
            Zbc_K(i,j).U3 = Zbc_K(i,j).U3(:,1:rank_min_K).';
        end
    end

end
