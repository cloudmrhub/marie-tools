function [x,r,k,resvec] = is_iter_gmres_vie(A,x,r,m,tol)

    [atype,afun,afcnstr] = iterchk(A);
    V = zeros(size(r,1),m+1,'like',x);
    V(:,1) = r / norm(r);

    H = zeros(m+1,m,'like',x);
    resvec = zeros(m,1,'like',x);

    for k = 1:m

        w = V(:,k);
        w = iterapp('mtimes',afun,atype,afcnstr,w);

        H(1:k,k) = V(:,1:k)' * w;
        w = w - V(:,1:k)*H(1:k,k);

        H(k+1,k) = norm(w);
        V(:,k+1) = w / H(k+1,k);

        rhs = zeros(k+1,1,'like',x);
        rhs(1) = norm(r);

        Hc = H(1:k+1,1:k);
        y = Hc \ rhs;
        res = rhs - Hc * y;
        resvec(k) = norm(res);

        if resvec(k) < tol
            break;
        end

    end

    resvec = resvec(1:k);
    xupd   = V(:,1:k)*y;
    x      = x + xupd;
    r      = V(:,1:k+1) * res;

end