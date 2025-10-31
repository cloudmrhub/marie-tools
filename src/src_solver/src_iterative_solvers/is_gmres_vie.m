function [x,iter,resvec] = is_gmres_vie(A,b,restart,tol,maxit,x0)

    [atype,afun,afcnstr] = iterchk(A);
    p = 1;
    if (nnz(x0))
        r = b - iterapp('mtimes',afun,atype,afcnstr,x0);
    else
        r = b;
    end

    % Calculate rhs norm
    bnorm = norm(b);
    clear b;

    resvec = zeros(restart*maxit,1,'like',r);
    it = 1;
    outit = 0;
    resvec(it) = norm(r)/bnorm;
    x = x0;

    while(resvec(it) > tol) && (outit < maxit)
        outit = outit+1;
        [x,r,p,resvec_inner] = is_iter_gmres_vie(A,x,r,restart,tol*bnorm);
        resvec(it+1:it+p) = resvec_inner/bnorm;
        it = it + p;
    end
    
    resvec = gather(resvec(1:it));
    iter = gather([p outit]);
    x = gather(x);
    clear r;

end