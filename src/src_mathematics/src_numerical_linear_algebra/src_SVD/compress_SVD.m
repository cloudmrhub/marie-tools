function[U,S,V] = compress_SVD(U,S,V,tol)
    s_norm  = diag(S)/S(1,1);
    sh_norm = s_norm(s_norm>tol);
    N       = length(sh_norm) + 1;
    
    if N > length(sh_norm)
        N = length(sh_norm);
    end

    S   = S(1:N,1:N);
    V   = V(:,1:N);
    U   = U(:,1:N);
end