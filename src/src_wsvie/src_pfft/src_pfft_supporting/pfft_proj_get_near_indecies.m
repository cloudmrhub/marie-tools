% function [idx_near_list] = pfft_proj_get_near_indecies(idx_x, idx_y, idx_z, N_near_1D, ext_n1, ext_n2, ext_n3)
% 
%     idx_x_min  = idx_x - ceil((N_near_1D-1)/2);
%     idx_y_min  = idx_y - ceil((N_near_1D-1)/2);
%     idx_z_min  = idx_z - ceil((N_near_1D-1)/2);
%     idx_x_max  = idx_x + ceil((N_near_1D-1)/2);
%     idx_y_max  = idx_y + ceil((N_near_1D-1)/2);
%     idx_z_max  = idx_z + ceil((N_near_1D-1)/2);
%     idx_cell_x = idx_x_min:idx_x_max;
%     idx_cell_y = idx_y_min:idx_y_max;
%     idx_cell_z = idx_z_min:idx_z_max;
%     idx_cell_x(idx_cell_x > ext_n1) = NaN;
%     idx_cell_y(idx_cell_y > ext_n2) = NaN;
%     idx_cell_z(idx_cell_z > ext_n3) = NaN;
%     idx_cell_x(idx_cell_x < 1) = NaN;
%     idx_cell_y(idx_cell_y < 1) = NaN;
%     idx_cell_z(idx_cell_z < 1) = NaN;
% 
%     [X,Y,Z] = meshgrid(idx_cell_x, idx_cell_y, idx_cell_z);
%     glob = sub2ind([ext_n1 ext_n2 ext_n3], X, Y, Z);
%     glob = permute(glob, [2 1 3]);
%     glob(isnan(glob)) = 0;
%     idx_near_list = glob(:).';
% 
% end
function idx_near_list = pfft_proj_get_near_indecies(idx_x, idx_y, idx_z, N_near_1D, ext_n1, ext_n2, ext_n3)

    r = ceil((N_near_1D-1)/2);
    idx_cell_x = (idx_x - r) : (idx_x + r);
    idx_cell_y = (idx_y - r) : (idx_y + r);
    idx_cell_z = (idx_z - r) : (idx_z + r);

    % Build grids with x as fastest-changing dim when linearized
    [X, Y, Z] = ndgrid(idx_cell_x, idx_cell_y, idx_cell_z);

    % Validity mask (no NaNs, within bounds)
    valid = (X >= 1 & X <= ext_n1) & ...
            (Y >= 1 & Y <= ext_n2) & ...
            (Z >= 1 & Z <= ext_n3);

    % allocate output grid of linear indices, zero where invalid
    glob = zeros(size(X), 'like', X);

    % only call sub2ind on valid subscripts
    if any(valid(:))
        glob(valid) = sub2ind([ext_n1 ext_n2 ext_n3], X(valid), Y(valid), Z(valid));
    end

    % row vector, zeros mark "no index"
    idx_near_list = glob(:).';
end
