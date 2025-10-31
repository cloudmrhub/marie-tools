function X = fast_pinv(A,rhs)

    [U,s,V] = svd(A,'econ');
    s = diag(s);
    s = 1./s(:);
    X = (V.*s.')*(transpose(conj(U(1,1:end))*rhs));
    
end
