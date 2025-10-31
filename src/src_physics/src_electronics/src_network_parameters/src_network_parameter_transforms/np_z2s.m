function s_params = np_z2s(z_params,Z0)
    % Transforms z parameters to s
    I = eye(size(z_params, 1));
    s_params = (z_params + Z0 * I) \ (z_params - Z0 * I);
end