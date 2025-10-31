function [ZR] = assembly_surf_ns_row(Ne_shield, index_shield, etod_shield, atod, L_obs, ro_x, ro_y, ro_z, rsX, rsY, rsZ, Np, W_t, Zi_j, A_o, F_o, first_node, second_node, k0, ZR, idx_coil, je_idx, L_j, L_src_vecs_j1, L_src_vecs_j2, j_1, j_2, ssj, Row_of_Interest)

    c = 1;
    store_ie = zeros(1,2);
    store_i = zeros(1,2);
    for ie = 1:Ne_shield
        for i=1:3
            if index_shield(atod(i,ie)) == Row_of_Interest
                store_ie(1,c) = ie;
                store_i(1,c) = i;
                c = c+1;
            end
        end
    end

    % Green's function and weight:
    roX    = reshape(ro_x(:,store_ie), [Np, 1, 2, 1]);
    roY    = reshape(ro_y(:,store_ie), [Np, 1, 2, 1]);
    roZ    = reshape(ro_z(:,store_ie), [Np, 1, 2, 1]);
    dx     = roX - rsX;
    dy     = roY - rsY;
    dz     = roZ - rsZ;
    Rmn_sq = dx.^2 + dy.^2 + dz.^2;
    Rmn    = sqrt(Rmn_sq);
    W_t_4D = repmat(W_t, 1, 1, size(Rmn,3), size(Rmn,4));
    W_GR   = W_t_4D .* exp(-1i * k0 * Rmn) ./ Rmn;

    for counter_row = 1:2
        %Green
        W_GR_in = squeeze(W_GR(:,:,counter_row,je_idx));
        % Compute for observer
        ie    = store_ie(1,counter_row);
        i     = store_i(1,counter_row);
        l_obs = squeeze(L_obs(:,:,ie));
        i_1   = second_node(i);
        i_2   = first_node(i);
        soi   = sign(etod_shield(i, ie));
        % Compute the lengths of the corresponding edges:
        L_i   = sqrt(l_obs(1,i)^2 + l_obs(2,i)^2 + l_obs(3,i)^2);
        % Compute dot products between the associated edge vectors:
        Li1_j1 = l_obs(:,i_1).' * L_src_vecs_j1; 
        Li1_j2 = l_obs(:,i_1).' * L_src_vecs_j2;
        Li2_j1 = l_obs(:,i_2).' * L_src_vecs_j1;
        Li2_j2 = l_obs(:,i_2).' * L_src_vecs_j2;
        % Retrieve the precomputed Zi_j integrals.
        Zi1_j1 = reshape(squeeze(Zi_j(:,:,i_1,j_1)),Np^2,[]);
        Zi1_j2 = reshape(squeeze(Zi_j(:,:,i_1,j_2)),Np^2,[]);
        Zi2_j1 = reshape(squeeze(Zi_j(:,:,i_2,j_1)),Np^2,[]);
        Zi2_j2 = reshape(squeeze(Zi_j(:,:,i_2,j_2)),Np^2,[]);

        % Assemble the Row
        ZZ          = Zi2_j2.*Li1_j1-Zi2_j1.*Li1_j2-Zi1_j2.*Li2_j1+Zi1_j1.*Li2_j2;
        ZZ          = reshape(ZZ, Np, Np, []);
        ZZf         = A_o.*ZZ+F_o;
        ZR_ie_je_ij = squeeze(soi*L_i*ssj.*L_j).*squeeze(sum(sum(W_GR_in.*ZZf)));
        ZR          = ZR + accumarray(idx_coil, ZR_ie_je_ij, [size(ZR,2), 1]).';

    end

end