function [idx_cell] = pfft_proj_find_expansion_cell(idx_x, idx_y, idx_z, N_exp_1D, n1, n2, n3)

    idx_range = -(N_exp_1D - 1) / 2:1:(N_exp_1D-1) / 2;

    idx_cell_x = idx_range(:) + idx_x;
    idx_cell_y = idx_range(:) + idx_y;
    idx_cell_z = idx_range(:) + idx_z;

    [X,Y,Z] = meshgrid(idx_cell_x, idx_cell_y, idx_cell_z);

    exp_glob_index = sub2ind([n1 n2 n3], X, Y, Z);

    exp_glob_index = permute(exp_glob_index, [2 1 3]);

    idx_cell = exp_glob_index(:).';

end