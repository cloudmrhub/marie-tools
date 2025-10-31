function [ZR] = assembly_surf_ns(index_shield, etod_shield, node_shield, elem_shield, index_coil, etod_coil, node_coil, elem_coil, ZR, Np_2D, k0)
    
    % Constants and quadrature
    first_node         = [3 1 2];
    second_node        = [2 3 1];
    [Np,wt,Z,z1,z2,z3] = gauss_2d_sie(Np_2D);
    W_t                = wt'*wt;
    Zi_j               = zeros(Np,Np,3,3);
    Rmn_x              = zeros(Np,Np);
    Rmn_y              = zeros(Np,Np);
    Rmn_z              = zeros(Np,Np);
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

    for ie = 1:Ne_shield
        ao    = squeeze(atod(:,ie));
        l_obs = squeeze(L_obs(:,:,ie));
    
        for je = 1:Ne_coil
            as    = squeeze(atsd(:,je));
            l_src = squeeze(L_src(:,:,je)); 
    
            % Green's function and weight:
            for ii=1:Np
                Rmn_x(ii,:) = ro_x(ii,ie)-rs_x(:,je).';
                Rmn_y(ii,:) = ro_y(ii,ie)-rs_y(:,je).';
                Rmn_z(ii,:) = ro_z(ii,ie)-rs_z(:,je).';
            end
            Rmn  = sqrt(Rmn_x.^2 + Rmn_y.^2 + Rmn_z.^2);
            W_GR = W_t .* exp(-1i*k0*Rmn)./Rmn;
            
            % -------------------------------------------
            % Loop over local edges in the shield and coil.
            % -------------------------------------------
            for i=1:3     % local edge index for shield (testing)
                for j=1:3 % local edge index for coil (basis)
                    % Get global indices for the shield and coil functions.
                    idx_shield = index_shield(ao(i));
                    idx_coil   = index_coil(as(j));
                    if (idx_shield ~= 0) && (idx_coil ~= 0)
                        % Orientation factors from etod:
                        soi = sign(etod_shield(i, ie));
                        ssj = sign(etod_coil(j, je));
                        % Determine alternative ordering using RWG convention:
                        i_1 = second_node(i);
                        i_2 = first_node(i);
                        j_1 = second_node(j);
                        j_2 = first_node(j);
                        % Compute dot products between the associated edge vectors:
                        Li1_j1 = l_obs(1,i_1)*l_src(1,j_1) + l_obs(2,i_1)*l_src(2,j_1) + l_obs(3,i_1)*l_src(3,j_1);
                        Li1_j2 = l_obs(1,i_1)*l_src(1,j_2) + l_obs(2,i_1)*l_src(2,j_2) + l_obs(3,i_1)*l_src(3,j_2);
                        Li2_j1 = l_obs(1,i_2)*l_src(1,j_1) + l_obs(2,i_2)*l_src(2,j_1) + l_obs(3,i_2)*l_src(3,j_1);
                        Li2_j2 = l_obs(1,i_2)*l_src(1,j_2) + l_obs(2,i_2)*l_src(2,j_2) + l_obs(3,i_2)*l_src(3,j_2);
                        % Compute the lengths of the corresponding edges:
                        L_i = sqrt(l_obs(1,i)^2 + l_obs(2,i)^2 + l_obs(3,i)^2);
                        L_j = sqrt(l_src(1,j)^2 + l_src(2,j)^2 + l_src(3,j)^2);
                        % Retrieve the precomputed Zi_j integrals.
                        Zi1_j1 = Zi_j(:,:,i_1,j_1);
                        Zi1_j2 = Zi_j(:,:,i_1,j_2);
                        Zi2_j1 = Zi_j(:,:,i_2,j_1);
                        Zi2_j2 = Zi_j(:,:,i_2,j_2);
                        % Combine the integrals in an RWG style:
                        ZZ = Li1_j1*Zi2_j2-Li1_j2*Zi2_j1-Li2_j1*Zi1_j2+Li2_j2*Zi1_j1;
                        % Weight with the constant terms:
                        ZZf = A_o.*ZZ+F_o;
                        % Multiply by the edge lengths and orientations, and integrate.
                        ZR_ie_je_ij = soi*ssj*L_i*L_j*(sum(sum(W_GR.*ZZf)));
                        % Update the global impedance matrix using the appropriate indices.
                        ZR(idx_shield, idx_coil) = ZR(idx_shield, idx_coil) + ZR_ie_je_ij;
                    end 
                end
            end
        end
    end
end