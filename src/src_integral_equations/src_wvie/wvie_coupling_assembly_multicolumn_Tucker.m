function[Zbw_N,Zbw_K,rank_min_N,rank_min_K] = wvie_coupling_assembly_multicolumn_Tucker(MREDM)
    
    k0   = MREDM.emc.k0;

    tol            = MREDM.inputs.tol_HOSVD_couple;
    Quad_order_vie = MREDM.inputs.Quad_order_vie;
    Quad_order_wie = MREDM.inputs.Quad_order_wie_coup;
    [wp_vie,z_vie] = gauss_1d(Quad_order_vie);
    [wp_wie,z_wie] = gauss_1d(Quad_order_wie);
    VIE_quads      = [Quad_order_vie; wp_vie; z_vie];
    WIE_quads      = [Quad_order_wie; wp_wie; z_wie];

    dims   = MREDM.dimensions;
    ql     = dims.ql;
    xd     = dims.r(:,:,:,1);
    yd     = dims.r(:,:,:,2);
    zd     = dims.r(:,:,:,3);
    n1     = dims.n1;
    n2     = dims.n2;
    n3     = dims.n3;
    res    = dims.res;
    N_wie  = dims.N_wie;
    Scoord = [xd(:) yd(:) zd(:)].';
    
    w_length = length(MREDM.WIE.coil.T_point);
    rp       = MREDM.WIE.coil.F_point(1:w_length,:).';
    r2       = MREDM.WIE.coil.S_point(1:w_length,:).';
    rn       = MREDM.WIE.coil.T_point(1:w_length,:).';

    Zbw_N = struct('G',repmat({[]},N_wie,ql), 'U1',repmat({[]},N_wie,ql), 'U2',repmat({[]},N_wie,ql), 'U3',repmat({[]},N_wie,ql));
    Zbw_K = struct('G',repmat({[]},N_wie,ql), 'U1',repmat({[]},N_wie,ql), 'U2',repmat({[]},N_wie,ql), 'U3',repmat({[]},N_wie,ql));

    if ql == 12
        for i = 1:N_wie
            [Zbw_N(i,1).G,Zbw_N(i,1).U1,Zbw_N(i,1).U2,Zbw_N(i,1).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_x(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_N(i,2).G,Zbw_N(i,2).U1,Zbw_N(i,2).U2,Zbw_N(i,2).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_x1(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,3).G,Zbw_N(i,3).U1,Zbw_N(i,3).U2,Zbw_N(i,3).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_x2(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,4).G,Zbw_N(i,4).U1,Zbw_N(i,4).U2,Zbw_N(i,4).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_x3(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,5).G,Zbw_N(i,5).U1,Zbw_N(i,5).U2,Zbw_N(i,5).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_y(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_N(i,6).G,Zbw_N(i,6).U1,Zbw_N(i,6).U2,Zbw_N(i,6).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_y1(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,7).G,Zbw_N(i,7).U1,Zbw_N(i,7).U2,Zbw_N(i,7).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_y2(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,8).G,Zbw_N(i,8).U1,Zbw_N(i,8).U2,Zbw_N(i,8).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_y3(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,9).G,Zbw_N(i,9).U1,Zbw_N(i,9).U2,Zbw_N(i,9).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_z(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_N(i,10).G,Zbw_N(i,10).U1,Zbw_N(i,10).U2,Zbw_N(i,10).U3] = hosvd(reshape(Assemble_tri_coupling_matrix_N_z1(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,11).G,Zbw_N(i,11).U1,Zbw_N(i,11).U2,Zbw_N(i,11).U3] = hosvd(reshape(Assemble_tri_coupling_matrix_N_z2(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_N(i,12).G,Zbw_N(i,12).U1,Zbw_N(i,12).U2,Zbw_N(i,12).U3] = hosvd(reshape(Assemble_tri_coupling_matrix_N_z3(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            
            [Zbw_K(i,1).G,Zbw_K(i,1).U1,Zbw_K(i,1).U2,Zbw_K(i,1).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_x(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_K(i,2).G,Zbw_K(i,2).U1,Zbw_K(i,2).U2,Zbw_K(i,2).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_x1(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,3).G,Zbw_K(i,3).U1,Zbw_K(i,3).U2,Zbw_K(i,3).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_x2(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,4).G,Zbw_K(i,4).U1,Zbw_K(i,4).U2,Zbw_K(i,4).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_x3(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,5).G,Zbw_K(i,5).U1,Zbw_K(i,5).U2,Zbw_K(i,5).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_y(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_K(i,6).G,Zbw_K(i,6).U1,Zbw_K(i,6).U2,Zbw_K(i,6).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_y1(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,7).G,Zbw_K(i,7).U1,Zbw_K(i,7).U2,Zbw_K(i,7).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_y2(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,8).G,Zbw_K(i,8).U1,Zbw_K(i,8).U2,Zbw_K(i,8).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_y3(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,9).G,Zbw_K(i,9).U1,Zbw_K(i,9).U2,Zbw_K(i,9).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_z(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_K(i,10).G,Zbw_K(i,10).U1,Zbw_K(i,10).U2,Zbw_K(i,10).U3] = hosvd(reshape(Assemble_tri_coupling_matrix_K_z1(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,11).G,Zbw_K(i,11).U1,Zbw_K(i,11).U2,Zbw_K(i,11).U3] = hosvd(reshape(Assemble_tri_coupling_matrix_K_z2(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
            [Zbw_K(i,12).G,Zbw_K(i,12).U1,Zbw_K(i,12).U2,Zbw_K(i,12).U3] = hosvd(reshape(Assemble_tri_coupling_matrix_K_z3(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0) * res^3,n1,n2,n3),tol);
        end
    elseif ql == 3
        for i = 1:N_wie
            [Zbw_N(i,1).G,Zbw_N(i,1).U1,Zbw_N(i,1).U2,Zbw_N(i,1).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_x(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_N(i,2).G,Zbw_N(i,2).U1,Zbw_N(i,2).U2,Zbw_N(i,2).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_y(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_N(i,3).G,Zbw_N(i,3).U1,Zbw_N(i,3).U2,Zbw_N(i,3).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_N_z(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_K(i,1).G,Zbw_K(i,1).U1,Zbw_K(i,1).U2,Zbw_K(i,1).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_x(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_K(i,2).G,Zbw_K(i,2).U1,Zbw_K(i,2).U2,Zbw_K(i,2).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_y(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
            [Zbw_K(i,3).G,Zbw_K(i,3).U1,Zbw_K(i,3).U2,Zbw_K(i,3).U3]     = hosvd(reshape(Assemble_tri_coupling_matrix_K_z(Scoord,rp(:,i),rn(:,i),r2(:,i),WIE_quads,VIE_quads,res,k0)  * res^3,n1,n2,n3),tol);
        end
    end

    % Crop to desired ranks
    rank_min_N = max([n1 n2 n3]);
    rank_min_K = max([n1 n2 n3]);
    for i = 1:N_wie
        for j = 1:ql
            [r1,r2,r3] = size(Zbw_N(i,j).G);
            if min([r1 r2 r3]) < rank_min_N
                rank_min_N = min([r1 r2 r3]);
            end
            [r1,r2,r3] = size(Zbw_K(i,j).G);
            if min([r1 r2 r3]) < rank_min_K
                rank_min_K = min([r1 r2 r3]);
            end
        end
    end
    for i = 1:N_wie
        for j = 1:ql
            Zbw_N(i,j).G = reshape(Zbw_N(i,j).G(1:rank_min_N,1:rank_min_N,1:rank_min_N),rank_min_N,rank_min_N^2);
            Zbw_N(i,j).U1 = Zbw_N(i,j).U1(:,1:rank_min_N);
            Zbw_N(i,j).U2 = Zbw_N(i,j).U2(:,1:rank_min_N);
            Zbw_N(i,j).U3 = Zbw_N(i,j).U3(:,1:rank_min_N).';
            Zbw_K(i,j).G = reshape(Zbw_K(i,j).G(1:rank_min_K,1:rank_min_K,1:rank_min_K),rank_min_K,rank_min_K^2);
            Zbw_K(i,j).U1 = Zbw_K(i,j).U1(:,1:rank_min_K);
            Zbw_K(i,j).U2 = Zbw_K(i,j).U2(:,1:rank_min_K);
            Zbw_K(i,j).U3 = Zbw_K(i,j).U3(:,1:rank_min_K).';
        end
    end

end
