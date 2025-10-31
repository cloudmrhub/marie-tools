function [ZR] = assembly_wire_surf_ns_column(Ne_coil, index_coil, etod_coil, atsd, L_src, roX_left, roY_left, roZ_left, roX_right, roY_right, roZ_right, rs_x, rs_y, rs_z, Np_1D, Np, W_t, Zi_j, A_o, F_o, first_node, second_node, k0, ZR, Nwie, Dl, Dr, fn_left, fn_right, tlp_left, tlp_right, Column_of_Interest)
    
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
    rsX        = reshape(rs_x(:,store_je), [1, Np, 1, 2]);
    rsY        = reshape(rs_y(:,store_je), [1, Np, 1, 2]);
    rsZ        = reshape(rs_z(:,store_je), [1, Np, 1, 2]);
    dx         = roX_left - rsX;
    dy         = roY_left - rsY;
    dz         = roZ_left - rsZ;
    Rmn_sq     = dx.^2 + dy.^2 + dz.^2;
    Rmn        = sqrt(Rmn_sq);
    W_t_4D     = repmat(W_t, 1, 1, size(Rmn,3), size(Rmn,4));
    W_GR_left  = W_t_4D .* exp(-1i * k0 * Rmn) ./ Rmn;
    dx         = roX_right - rsX;
    dy         = roY_right - rsY;
    dz         = roZ_right - rsZ;
    Rmn_sq     = dx.^2 + dy.^2 + dz.^2;
    Rmn        = sqrt(Rmn_sq);
    W_GR_right = W_t_4D .* exp(-1i * k0 * Rmn) ./ Rmn;
    
    for counter_column = 1:2

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
        Lleft_j1  = l_src(:,j_1).' * ( tlp_left .* fn_left ); 
        Lleft_j2  = l_src(:,j_2).' * ( tlp_left .* fn_left );
        Lright_j1 = l_src(:,j_1).' * ( tlp_right.* fn_right); 
        Lright_j2 = l_src(:,j_2).' * ( tlp_right.* fn_right);

        % Retrieve the precomputed Zi_j integrals.
        Zi1_j1 = reshape(squeeze(Zi_j(:,:,ones(Nwie,1),j_1)),Np_1D*Np,[]);
        Zi1_j2 = reshape(squeeze(Zi_j(:,:,ones(Nwie,1),j_2)),Np_1D*Np,[]);

        % Assemble the Column
        ZZ_left     = -Zi1_j2.*Lleft_j2+Zi1_j1.*Lleft_j1;
        ZZ_right    = -Zi1_j2.*Lright_j2+Zi1_j1.*Lright_j1;
        ZZ_left     = reshape(ZZ_left, Np_1D, Np, []);
        ZZ_right    = reshape(ZZ_right, Np_1D, Np, []);
        ZZf_left    = A_o.*ZZ_left+F_o;
        ZZf_right   = A_o.*ZZ_right+F_o;
    
        ZR_left_je_j  = squeeze(ssj.*L_j).*Dl/2.*squeeze(sum(sum(squeeze(W_GR_left(:,:,1:Nwie,counter_column)).*ZZf_left)));
        ZR_right_je_j = squeeze(ssj.*L_j).*Dr/2.*squeeze(sum(sum(squeeze(W_GR_right(:,:,1:Nwie,counter_column)).*ZZf_right)));

        ZR = ZR + ZR_left_je_j + ZR_right_je_j;
    end

end