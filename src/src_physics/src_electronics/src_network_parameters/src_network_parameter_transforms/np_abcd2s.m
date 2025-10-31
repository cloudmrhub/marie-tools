function s_params = np_abcd2s(abcd_params,z0)
    % Transforms abcd parameters to s
    siz = size(abcd_params);
    m = siz(1)/2; 
    if siz(1) == 2
        [a, b, c, d] = deal(abcd_params(1,1), abcd_params(1,2) ./ z0, abcd_params(2,1) .* z0, abcd_params(2,2));
        delta = a+b+c+d;
        s_params = [(a+b-c-d), 2*(a.*d-b.*c); 2*ones(size(b)), (-a+b-c+d)] ./ repmat(delta, [2 2 1]);
    else
        I = eye(m);
        A = abcd_params(1:m,1:m);
        B = abcd_params(1:m,(m+1):(2*m));
        C = abcd_params((m+1):(2*m),1:m);
        D = abcd_params((m+1):(2*m),(m+1):(2*m));
        BB = B/(z0*I);
        CC = C*(z0*I);
        denom = A + BB + CC + D;
        S12 = ((A-BB-CC+D)-(A+BB-CC-D)/(A+BB+CC+D)*(A-BB+CC-D))/2;
        s_params = [(A + BB - CC - D)/denom, S12; 2*I/denom, denom\(-A + BB - CC + D)];
    end

end
