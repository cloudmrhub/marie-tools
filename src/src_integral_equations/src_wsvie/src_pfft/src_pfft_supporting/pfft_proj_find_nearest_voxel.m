function [i_vie, j_vie, k_vie] = pfft_proj_find_nearest_voxel(r_sie, xd, yd, zd, res)

    [list_idx_x] = find((xd >= (r_sie(1) - res)) & (xd <= (r_sie(1) + res)));
    [list_idx_y] = find((yd >= (r_sie(2) - res)) & (yd <= (r_sie(2) + res)));
    [list_idx_z] = find((zd >= (r_sie(3) - res)) & (zd <= (r_sie(3) + res)));

    for i_x = 1:length(list_idx_x)
        for j_y = 1:length(list_idx_y)
            for k_z = 1:length(list_idx_z)
                dist = norm(r_sie - [xd(list_idx_x(i_x)); yd(list_idx_y(j_y)); zd(list_idx_z(k_z))],Inf);
                if dist - res/2 < eps
                    i_vie = list_idx_x(i_x);
                    j_vie = list_idx_y(j_y);
                    k_vie = list_idx_z(k_z);
                    return;
                end
            end
        end
    end
    
end

