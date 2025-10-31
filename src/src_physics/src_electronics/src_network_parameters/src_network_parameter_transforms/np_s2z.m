function z_params = np_s2z(s_params,Z0)
    % Transforms s parameters to z
    I = eye(size(s_params, 1));
    z_params = (Z0 * (I + s_params)) /(I - s_params);
end