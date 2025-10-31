function [E] =  pfft_proj_pwx_to_collocation(vox_center, res, col_points, GL_order, pwx_sz, N_basis, k0, ce) 
    
    [wt_vie,z]  = gauss_1d(GL_order); 
    [rows,cols] = size(col_points);
    N           = size(vox_center,2);
    M           = rows * cols;
    E           = zeros(M,pwx_sz*N);
    % For all projection voxels
    parfor voxel = 1:N
        % For all collocation points
        for collocation_point = 1:cols
            r_point = col_points(:,collocation_point);
            % For all basis functions components
            for basis = 1:N_basis
                E_cur = zeros(3, 3); % Green function kernel
                for i1_vie = 1:length(wt_vie)
                    % For every Gaussian point per x
                    x_src = vox_center(1, voxel) + res / 2 * z(i1_vie);
                    for j2_vie = 1:length(wt_vie)
                        % For every Gaussian point per y
                        y_src = vox_center(2, voxel) + res / 2 * z(j2_vie);
                        for k3_vie = 1:length(wt_vie)
                            % For every Gaussian point per z
                            z_src = vox_center(3, voxel) + res / 2 * z(k3_vie);
                            r_src = [x_src; y_src; z_src]; % Coordinates
                            Weight = res^3 * wt_vie(i1_vie) * wt_vie(j2_vie) * wt_vie(k3_vie) / 8.0; % Weights
                            R_vec =  r_src - r_point; % Distance vector
                            R     =  norm(R_vec); % Distance
                            green = exp(-1i * k0 * R) / (4 * pi * R^3);
                            P     = (3 + 3 * 1i * k0 * R - k0^2 * R^2) / R^2;
                            Q     = (1 + 1i * k0 * R - k0^2 * R^2) / P;
                            Gxx   = R_vec(1)^2 - Q;
                            Gxy   = R_vec(1)*R_vec(2);
                            Gxz   = R_vec(1)*R_vec(3);
                            Gyy   = R_vec(2)^2 - Q;
                            Gyz   = R_vec(2)*R_vec(3);
                            Gzz   = R_vec(3)^2 - Q;
                            DGF = [Gxx Gxy Gxz;Gxy Gyy Gyz;Gxz Gyz Gzz];
                            switch basis
                                case 1
                                    E_cur(:,:) = E_cur(:,:) + Weight.*(green*P).*DGF;
                                case 2
                                    E_cur(:,:) = E_cur(:,:) + Weight.*(green*P).*DGF.* z(i1_vie) / 2;
                                case 3
                                    E_cur(:,:) = E_cur(:,:) + Weight.*(green*P).*DGF.* z(j2_vie) / 2;
                                case 4
                                    E_cur(:,:) = E_cur(:,:) + Weight.*(green*P).*DGF.* z(k3_vie) / 2;
                            end
                        end
                    end
                end
                i_x             = collocation_point;
                i_y             = collocation_point + cols;
                i_z             = collocation_point + 2 * cols;
                j_x             = voxel + (basis-1) * N;
                j_y             = j_x + N * N_basis;
                j_z             = j_y + N * N_basis;
                E_temp          = zeros(M,pwx_sz * N);
                E_temp(i_x,j_x) = E_cur(1,1);
                E_temp(i_y,j_x) = E_cur(2,1);
                E_temp(i_z,j_x) = E_cur(3,1);
                E_temp(i_x,j_y) = E_cur(1,2);
                E_temp(i_y,j_y) = E_cur(2,2);
                E_temp(i_z,j_y) = E_cur(3,2);
                E_temp(i_x,j_z) = E_cur(1,3);
                E_temp(i_y,j_z) = E_cur(2,3);
                E_temp(i_z,j_z) = E_cur(3,3);
                E               = E + E_temp;
            end
        end
    end
    E = 1 / ce * E;
end
