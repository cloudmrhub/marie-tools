function [U,V] = assembly_surf_ns_aca_fun(index_shield, etod_shield, node_shield, elem_shield, index_coil, etod_coil, node_coil, elem_coil, Np_2D, k0, tol)
    
    % Constants and quadrature
    first_node         = [3 1 2];
    second_node        = [2 3 1];
    [Np,wt,Z,z1,z2,z3] = gauss_2d_sie(Np_2D);
    W_t                = wt'*wt;
    Zi_j               = zeros(Np,Np,3,3);
    A_o                = 1i*k0*ones(Np,Np);
    F_o                = 1/(1i*k0)*4*ones(Np,Np);
    
    % Precompute the Zi_j integrals (for all combinations of quadrature points)
    % Dimensions: Np x Np x 3 x 3.
    for ii=1:3
        for jj=1:3
            Zi_j(:,:,ii,jj) = Z(:,ii)*Z(:,jj)';
        end
    end

    % Enumerate all NS interactions between shield and coil:
    % For shield–coil interactions, we consider every shield element with every
    % coil element.    Ne_coil   = size(elem_coil, 2);
    Ne_coil   = size(elem_coil, 2);
    Ne_shield = size(elem_shield, 2);

    % Entries for shield points
    atod         = abs(etod_shield);     
    node_test_1  = elem_shield(1,:);
    node_test_2  = elem_shield(2,:);
    node_test_3  = elem_shield(3,:);
    ro_1         = node_shield(:,node_test_1);
    ro_2         = node_shield(:,node_test_2);
    ro_3         = node_shield(:,node_test_3);
    L_obs        = zeros(3,3,Ne_shield);
    L_obs(:,1,:) = ro_2-ro_3;
    L_obs(:,2,:) = ro_3-ro_1;
    L_obs(:,3,:) = ro_1-ro_2;
    ro_x         = z1.'*ro_1(1,:)+z2.'*ro_2(1,:)+z3.'*ro_3(1,:);
    ro_y         = z1.'*ro_1(2,:)+z2.'*ro_2(2,:)+z3.'*ro_3(2,:);
    ro_z         = z1.'*ro_1(3,:)+z2.'*ro_2(3,:)+z3.'*ro_3(3,:);
    roX          = reshape(ro_x, [Np, 1, Ne_shield, 1]); 
    roY          = reshape(ro_y, [Np, 1, Ne_shield, 1]);
    roZ          = reshape(ro_z, [Np, 1, Ne_shield, 1]);

    % Entries for coil points
    atsd         = abs(etod_coil);
    node_basis_1 = elem_coil(1,:);
    node_basis_2 = elem_coil(2,:);
    node_basis_3 = elem_coil(3,:);
    rs_1         = node_coil(:,node_basis_1);
    rs_2         = node_coil(:,node_basis_2);
    rs_3         = node_coil(:,node_basis_3);
    L_src        = zeros(3,3,Ne_coil);
    L_src(:,1,:) = rs_2-rs_3;
    L_src(:,2,:) = rs_3-rs_1;
    L_src(:,3,:) = rs_1-rs_2;
    rs_x         = z1.'*rs_1(1,:)+z2.'*rs_2(1,:)+z3.'*rs_3(1,:);
    rs_y         = z1.'*rs_1(2,:)+z2.'*rs_2(2,:)+z3.'*rs_3(2,:);
    rs_z         = z1.'*rs_1(3,:)+z2.'*rs_2(3,:)+z3.'*rs_3(3,:);
    rsX          = reshape(rs_x, [1, Np, 1, Ne_coil]); 
    rsY          = reshape(rs_y, [1, Np, 1, Ne_coil]);
    rsZ          = reshape(rs_z, [1, Np, 1, Ne_coil]);

    % Precompute for Rows
    idx_coil = [];
    j_coil   = [];
    ssj      = [];
    j_1      = [];
    j_2      = [];
    je_idx   = [];
    for je = 1:Ne_coil
        as = squeeze(atsd(:,je));
        for j=1:3
            if index_coil(as(j)) ~= 0
                j_coil   = [j_coil;j];
                idx_coil = [idx_coil;index_coil(as(j))];
                ssj      = [ssj;sign(etod_coil(j, je))];
                j_1      = [j_1;second_node(j)];
                j_2      = [j_2;first_node(j)];
                je_idx   = [je_idx; je];
            end
        end
    end
    global_idx1   = (je_idx - 1) * 3 + j_coil; 
    global_idx_j1 = (je_idx - 1) * 3 + j_1;
    global_idx_j2 = (je_idx - 1) * 3 + j_2;
    L_src_flat    = reshape(L_src, 3, []); 
    L_src_vecs    = L_src_flat(:, global_idx1);  
    L_j           = sqrt(sum(L_src_vecs.^2, 1)).'; 
    L_src_vecs_j1 = L_src_flat(:, global_idx_j1); 
    L_src_vecs_j2 = L_src_flat(:, global_idx_j2); 

    % Precompute for Columns
    idx_shield = [];
    i_shield   = [];
    soi        = [];
    i_1        = [];
    i_2        = [];
    ie_idx     = [];
    for ie = 1:Ne_shield
        ao = squeeze(atod(:,ie));
        for i=1:3
            if index_shield(ao(i)) ~= 0
                i_shield   = [i_shield;i];
                idx_shield = [idx_shield;index_shield(ao(i))];
                soi      = [soi;sign(etod_shield(i, ie))];
                i_1      = [i_1;second_node(i)];
                i_2      = [i_2;first_node(i)];
                ie_idx   = [ie_idx; ie];
            end
        end
    end
    global_idx2   = (ie_idx - 1) * 3 + i_shield; 
    global_idx_i1 = (ie_idx - 1) * 3 + i_1;
    global_idx_i2 = (ie_idx - 1) * 3 + i_2;
    L_obs_flat    = reshape(L_obs, 3, []); 
    L_obs_vecs    = L_obs_flat(:, global_idx2);  
    L_i           = sqrt(sum(L_obs_vecs.^2, 1)).'; 
    L_obs_vecs_i1 = L_obs_flat(:, global_idx_i1); 
    L_obs_vecs_i2 = L_obs_flat(:, global_idx_i2); 

    % Initialize ACA
    I = 1;
	U = [];
	V = [];

    ZR_row = zeros(1,max(index_coil));
    ZR_col = zeros(max(index_shield),1);
	
	%Ri = zeros(M, 1);
	SF2 = 0;

    for k=1:min(Ne_shield,Ne_coil)
        % Compute rows
        Ri = assembly_surf_ns_row(Ne_shield, index_shield, etod_shield, atod, L_obs, ro_x, ro_y, ro_z, rsX, rsY, rsZ, Np, W_t, Zi_j, A_o, F_o, first_node, second_node, k0, ZR_row, idx_coil, je_idx, L_j, L_src_vecs_j1, L_src_vecs_j2, j_1, j_2, ssj, I);
        if k > 1
            Ri = Ri - U(I,:)*V';
        end
        [~,J] = max(abs(Ri));
        Vk = (Ri / Ri(J))';
        
        % Compute columns
        Uk = assembly_surf_ns_column(Ne_coil, index_coil, etod_coil, atsd, L_src, roX, roY, roZ, rs_x, rs_y, rs_z, Np, W_t, Zi_j, A_o, F_o, first_node, second_node, k0, ZR_col, idx_shield, ie_idx, L_i, L_obs_vecs_i1, L_obs_vecs_i2, i_1, i_2, soi, J);
        if k > 1
            Uk = Uk - U*(V(J,:))';
        end

        % Compute Approximation
        nVk = norm(Vk);
        nUk = norm(Uk);
        SF2 = SF2 + (nVk*nUk)^2; 
        if k > 1
            SF2 = SF2 + 2*sum(real((U'*Uk).*(V'*Vk))); 
        end

        % Evaluate if we are done
        U = [U Uk];
        V = [V Vk];
        if nUk*nVk <= tol*sqrt(SF2)
            break;
        end
        Uk = abs(Uk);
        Uk(I) = 0;
        [~,I] = max(Uk);
        
    end
    
end
