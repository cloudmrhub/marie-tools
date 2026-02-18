function abcd_params = np_s2abcd(s_params, z0)
    % Transforms s parameters to abcd
    siz = size(s_params);
    m = siz(1)/2;
    if siz(1) == 2
        [s11, s12, s21, s22] = deal(s_params(1,1), s_params(1,2), s_params(2,1), s_params(2,2));
        temp1 = 1 - s11;        temp2 = 1 + s11;
        temp3 = 1 - s22;        temp4 = 1 + s22;
        s12s21 = s12 .* s21;
        abcd_params = [temp2 .* temp3 + s12s21, (temp2 .* temp4 - s12s21) .* z0; (temp1 .* temp3 - s12s21) ./ z0, temp1 .* temp4 + s12s21] ./ repmat(2 * s21, [2 2 1]);    
    else
        I = eye(m);
        s11 = s_params(1:m,1:m);
        s12 = s_params(1:m,(m+1):(2*m));
        s21 = s_params((m+1):(2*m),1:m);
        s22 = s_params((m+1):(2*m),(m+1):(2*m));
        A = 0.5*((I+s11)/s21*(I-s22)+s12);
        B = z0*0.5*((I+s11)/s21*(I+s22)-s12);
        C = 0.5/z0*((I-s11)/s21*(I-s22)-s12);
        D = 0.5*((I-s11)/s21*(I+s22)+s12);
        abcd_params = [A , B ; C , D];
    end
end