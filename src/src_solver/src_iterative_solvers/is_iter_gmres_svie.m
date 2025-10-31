function [x,r,k,resvec] = is_iter_gmres_svie(A,x,r,m,L1,pL_L,L2,Mbprec,N_sie,tol,rnorm)

    [atype,afun,afcnstr] = iterchk(A);
    V = zeros(size(r,1),m+1,'like',Mbprec);
    V(:,1) = r / rnorm;  

    H = zeros(m+1,m,'like',Mbprec);
    resvec = zeros(m,1,'like',Mbprec);

    for k = 1:m

        w = iterapp('mtimes',afun,atype,afcnstr,V(:,k));

        w1 = w(1:N_sie);
        w2 = w(N_sie+1:end);

        if isa(L1, 'function_handle')  
            w1  = L1(w1);
        elseif ( isa(L1,'gpuArray') && isa(Mbprec,'gpuArray') ) || ( ~isa(L1,'gpuArray') && ~isa(Mbprec,'gpuArray') )
            w1 = L2 * (L1 * w1(pL_L));
        else
            w1 = gpuArray(L2 * (L1 * gather(w1(pL_L))));
        end
        w2 = Mbprec .* w2;
        
        w  = [w1;w2];

        H(1:k,k) = V(:,1:k)' * w; 
        w = w - V(:,1:k)*H(1:k,k);

        H(k+1,k) = norm(w); 
        V(:,k+1) = w / H(k+1,k);

        rhs = zeros(k+1,1,'like',Mbprec);
        rhs(1) = rnorm;

        Hc = H(1:k+1,1:k);
        y = fast_pinv(Hc,rnorm);
        res = rhs - Hc * y;
        resvec(k) = norm(res);

        if resvec(k) < tol
            break;
        end

    end

    resvec = resvec(1:k);
    x      = x + V(:,1:k)*y;
    r      = V(:,1:k+1) * res;

    clear L1 L2 Mbprec;

end
