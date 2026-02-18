function y_params = np_s2y(s_params,Z0)
    % Transforms s parameters to y
    I = eye(size(s_params, 1));
    y_params = (I - s_params) / (Z0 * (I + s_params));
end
