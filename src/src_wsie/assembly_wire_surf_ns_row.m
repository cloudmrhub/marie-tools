function [ZR] = assembly_wire_surf_ns_row(Dl, Dr, fn_left, fn_right, tlp_left, tlp_right, ro_x_left, ro_y_left, ro_z_left, ro_x_right, ro_y_right, ro_z_right, rsX, rsY, rsZ, Np_1D, Np, W_t, Zi_j, A_o, F_o, k0, ZR, idx_coil, je_idx, L_j, L_src_vecs_j1, L_src_vecs_j2, j_1, j_2, ssj, Row_of_Interest)

    % Green's function and weight for left segment
    dx         = reshape(ro_x_left(:,Row_of_Interest), [Np_1D, 1, 1, 1]) - rsX;
    dy         = reshape(ro_y_left(:,Row_of_Interest), [Np_1D, 1, 1, 1]) - rsY;
    dz         = reshape(ro_z_left(:,Row_of_Interest), [Np_1D, 1, 1, 1]) - rsZ;
    Rmn_sq     = dx.^2 + dy.^2 + dz.^2;
    Rmn        = sqrt(Rmn_sq);
    W_t_4D     = repmat(W_t, 1, 1, size(Rmn,3), size(Rmn,4));
    W_GR_left  = W_t_4D .* exp(-1i * k0 * Rmn) ./ Rmn;
    W_GR_left  = squeeze(W_GR_left(:,:,:,je_idx));

    % Green's function and weight for right segment:
    dx         = reshape(ro_x_right(:,Row_of_Interest), [Np_1D, 1, 1, 1]) - rsX;
    dy         = reshape(ro_y_right(:,Row_of_Interest), [Np_1D, 1, 1, 1]) - rsY;
    dz         = reshape(ro_z_right(:,Row_of_Interest), [Np_1D, 1, 1, 1]) - rsZ;
    Rmn_sq     = dx.^2 + dy.^2 + dz.^2;
    Rmn        = sqrt(Rmn_sq);
    W_GR_right = W_t_4D .* exp(-1i * k0 * Rmn) ./ Rmn;
    W_GR_right = squeeze(W_GR_right(:,:,:,je_idx));

    % Compute dot products between the associated edge vectors:
    Lleft_j1  = ( fn_left(:,Row_of_Interest)  * tlp_left(:,Row_of_Interest).'  ) * L_src_vecs_j1; 
    Lleft_j2  = ( fn_left(:,Row_of_Interest)  * tlp_left(:,Row_of_Interest).'  ) * L_src_vecs_j2;
    Lright_j1 = ( fn_right(:,Row_of_Interest) * tlp_right(:,Row_of_Interest).' ) * L_src_vecs_j1; 
    Lright_j2 = ( fn_right(:,Row_of_Interest) * tlp_right(:,Row_of_Interest).' ) * L_src_vecs_j2;

    % Retrieve the precomputed Zi_j integrals.
    Zi1_j1 = reshape(squeeze(Zi_j(:,:,j_1)),Np_1D*Np,[]);
    Zi1_j2 = reshape(squeeze(Zi_j(:,:,j_2)),Np_1D*Np,[]);

    % Assemble the Row
    ZZ_left     = -Zi1_j2.*Lleft_j2+Zi1_j1.*Lleft_j1;
    ZZ_right    = -Zi1_j2.*Lright_j2+Zi1_j1.*Lright_j1;
    ZZ_left     = reshape(ZZ_left, Np_1D, Np, []);
    ZZ_right    = reshape(ZZ_right, Np_1D, Np, []);
    ZZf_left    = A_o.*ZZ_left+F_o;
    ZZf_right   = A_o.*ZZ_right+F_o;

    ZR_left_je_j  = Dl(Row_of_Interest)/2 * squeeze(ssj.*L_j).*squeeze(sum(sum(W_GR_left.*ZZf_left)));
    ZR_right_je_j = Dr(Row_of_Interest)/2 * squeeze(ssj.*L_j).*squeeze(sum(sum(W_GR_right.*ZZf_right)));

    ZR            = ZR + accumarray(idx_coil, ZR_left_je_j, [size(ZR,2), 1]).' ...
                       + accumarray(idx_coil, ZR_right_je_j, [size(ZR,2), 1]).';
end