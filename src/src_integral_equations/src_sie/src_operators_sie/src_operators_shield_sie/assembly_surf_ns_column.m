function [ZR] = assembly_surf_ns_column(Ne_coil, index_coil, etod_coil, atsd, L_src, roX, roY, roZ, rs_x, rs_y, rs_z, Np, W_t, Zi_j, A_o, F_o, first_node, second_node, k0, ZR, idx_shield, ie_idx, L_i, L_obs_vecs_i1, L_obs_vecs_i2, i_1, i_2, soi, Column_of_Interest)
    
    c = 1;
    store_je = zeros(1,2);
    store_j = zeros(1,2);
    for je = 1:Ne_coil
        for j=1:3
            if index_coil(atsd(j,je)) == Column_of_Interest
                store_je(1,c) = je;
                store_j(1,c) = j;
                c = c+1;
            end
        end
    end

    % Green's function and weight:
    rsX    = reshape(rs_x(:,store_je), [1, Np, 1, 2]);
    rsY    = reshape(rs_y(:,store_je), [1, Np, 1, 2]);
    rsZ    = reshape(rs_z(:,store_je), [1, Np, 1, 2]);
    dx     = roX - rsX;
    dy     = roY - rsY;
    dz     = roZ - rsZ;
    Rmn_sq = dx.^2 + dy.^2 + dz.^2;
    Rmn    = sqrt(Rmn_sq);
    W_t_4D = repmat(W_t, 1, 1, size(Rmn,3), size(Rmn,4));
    W_GR   = W_t_4D .* exp(-1i * k0 * Rmn) ./ Rmn;
    
    for counter_column = 1:2
        %Green
        W_GR_in = squeeze(W_GR(:,:,ie_idx,counter_column));
        % Compute for source
        je    = store_je(1,counter_column);
        j     = store_j(1,counter_column);
        l_src = squeeze(L_src(:,:,je));
        j_1   = second_node(j);
        j_2   = first_node(j);
        ssj   = sign(etod_coil(j, je));
        % Compute the lengths of the corresponding edges:
        L_j   = sqrt(l_src(1,j)^2 + l_src(2,j)^2 + l_src(3,j)^2);
        % Compute dot products between the associated edge vectors:
        Li1_j1 = l_src(:,j_1).' * L_obs_vecs_i1; 
        Li1_j2 = l_src(:,j_2).' * L_obs_vecs_i1;
        Li2_j1 = l_src(:,j_1).' * L_obs_vecs_i2;
        Li2_j2 = l_src(:,j_2).' * L_obs_vecs_i2;
        % Retrieve the precomputed Zi_j integrals.
        Zi1_j1 = reshape(squeeze(Zi_j(:,:,i_1,j_1)),Np^2,[]);
        Zi1_j2 = reshape(squeeze(Zi_j(:,:,i_1,j_2)),Np^2,[]);
        Zi2_j1 = reshape(squeeze(Zi_j(:,:,i_2,j_1)),Np^2,[]);
        Zi2_j2 = reshape(squeeze(Zi_j(:,:,i_2,j_2)),Np^2,[]);

        % Assemble the Column
        ZZ          = Zi2_j2.*Li1_j1-Zi2_j1.*Li1_j2-Zi1_j2.*Li2_j1+Zi1_j1.*Li2_j2;
        ZZ          = reshape(ZZ, Np, Np, []);
        ZZf         = A_o.*ZZ+F_o;
        ZR_ie_je_ij = squeeze(ssj*L_j*soi.*L_i).*squeeze(sum(sum(W_GR_in.*ZZf)));
        ZR          = ZR + accumarray(idx_shield, ZR_ie_je_ij, [size(ZR,1), 1]);
    end

end