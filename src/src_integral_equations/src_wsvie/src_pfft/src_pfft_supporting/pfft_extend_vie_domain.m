function [r_ext, ext_mask] = pfft_extend_vie_domain(x_3d, y_3d, z_3d, x, y, z, res, idxS, mask, xc, yc, zc, d_coil2bound, N_exp_half)

    %% Build bounding box of BODY
    Xb_bb   = zeros(3,2);
    Xb_bb(1,1) = min(x_3d(idxS)) - res;   Xb_bb(1,2) = max(x_3d(idxS)) + res;
    Xb_bb(2,1) = min(y_3d(idxS)) - res;   Xb_bb(2,2) = max(y_3d(idxS)) + res;
    Xb_bb(3,1) = min(z_3d(idxS)) - res;   Xb_bb(3,2) = max(z_3d(idxS)) + res;

    %% Build bounding box of COIL
    Xc_bb   = zeros(3,2);
    Xc_bb(1,1) = min(xc);   Xc_bb(1,2) = max(xc);
    Xc_bb(2,1) = min(yc);   Xc_bb(2,2) = max(yc);
    Xc_bb(3,1) = min(zc);   Xc_bb(3,2) = max(zc);

    %% Build coil-expanded bounding box Xext_bb
    Xext_bb = zeros(3,2);

    for i = 1:3
        % left side
        if Xb_bb(i,1) < Xc_bb(i,1) - d_coil2bound
            Xext_bb(i,1) = Xb_bb(i,1);
        else
            n_shift = round(abs(Xb_bb(i,1) - Xc_bb(i,1)) / res);
            Xext_bb(i,1) = Xb_bb(i,1) - res * (n_shift + N_exp_half);
        end

        % right side
        if Xb_bb(i,2) > Xc_bb(i,2) + d_coil2bound
            Xext_bb(i,2) = Xb_bb(i,2);
        else
            n_shift = round(abs(Xb_bb(i,2) - Xc_bb(i,2)) / res);
            Xext_bb(i,2) = Xb_bb(i,2) + res * (n_shift + N_exp_half);
        end
    end

    %% Merge extended bounding box with the original VIE grid extents
    xmin_body = x(1);  xmax_body = x(end);
    ymin_body = y(1);  ymax_body = y(end);
    zmin_body = z(1);  zmax_body = z(end);

    xmin_ext = Xext_bb(1,1);  xmax_ext = Xext_bb(1,2);
    ymin_ext = Xext_bb(2,1);  ymax_ext = Xext_bb(2,2);
    zmin_ext = Xext_bb(3,1);  zmax_ext = Xext_bb(3,2);

    % final bounds while preserving voxel alignment
    xmin = min(xmin_body, xmin_ext);
    xmax = max(xmax_body, xmax_ext);
    ymin = min(ymin_body, ymin_ext);
    ymax = max(ymax_body, ymax_ext);
    zmin = min(zmin_body, zmin_ext);
    zmax = max(zmax_body, zmax_ext);

    %% Snap boundaries relative to BODY grid (prevents 1-voxel offset)
    xmin = xmin_body + res * floor((xmin - xmin_body)/res);
    xmax = xmax_body + res * ceil ((xmax - xmax_body)/res);

    ymin = ymin_body + res * floor((ymin - ymin_body)/res);
    ymax = ymax_body + res * ceil ((ymax - ymax_body)/res);

    zmin = zmin_body + res * floor((zmin - zmin_body)/res);
    zmax = zmax_body + res * ceil ((zmax - zmax_body)/res);

    %% Build voxel-aligned extended grid
    Nx = round((xmax - xmin) / res);
    Ny = round((ymax - ymin) / res);
    Nz = round((zmax - zmin) / res);

    x_ext = xmin + (0:Nx)*res;
    y_ext = ymin + (0:Ny)*res;
    z_ext = zmin + (0:Nz)*res;

    [Xext, Yext, Zext] = ndgrid(x_ext, y_ext, z_ext);
    r_ext = cat(4, Xext, Yext, Zext);

    %% Map original mask into extended domain
    xb_min = min(x);  xb_max = max(x);
    yb_min = min(y);  yb_max = max(y);
    zb_min = min(z);  zb_max = max(z);

    idx_vie_x = (x_ext >= xb_min - res/3) & (x_ext <= xb_max + res/3);
    idx_vie_y = (y_ext >= yb_min - res/3) & (y_ext <= yb_max + res/3);
    idx_vie_z = (z_ext >= zb_min - res/3) & (z_ext <= zb_max + res/3);

    ext_mask = zeros(length(x_ext), length(y_ext), length(z_ext));
    ext_mask(idx_vie_x, idx_vie_y, idx_vie_z) = mask;
end
