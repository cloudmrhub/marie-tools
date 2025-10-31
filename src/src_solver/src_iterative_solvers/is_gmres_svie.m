function [x,iter,resvec] = is_gmres_svie(A,b,restart,tol,maxit,L1,pL_L,L2,Mbprec,N_sie,x0)

    [atype,afun,afcnstr] = iterchk(A);

    if (nnz(x0))
        r = b - iterapp('mtimes',afun,atype,afcnstr,x0);
    else
        r = b;
    end
    bp = b;

    r1  = r(1:N_sie);
    bp1 = bp(1:N_sie);
    r2  = r(N_sie+1:end);
    bp2 = bp(N_sie+1:end);

    if ( isa(L1,'gpuArray') && isa(Mbprec,'gpuArray') ) || ( ~isa(L1,'gpuArray') && ~isa(Mbprec,'gpuArray') )
        r1  = L2 * (L1 * r1(pL_L));
        bp1 = L2 * (L1 * bp1(pL_L));
    else
        r1  = gpuArray(L2 * (L1 * gather(r1(pL_L))));
        bp1 = gpuArray(L2 * (L1 * gather(bp1(pL_L))));
    end
    r2  = Mbprec .* r2;
    bp2 = Mbprec .* bp2;
    
    r   = [r1;r2];
    bp  = [bp1;bp2];

    bnorm = norm(bp);
    clear bp bp1 bp2 r1 r2;

    resvec = zeros(restart*maxit,1,'like',Mbprec);
    it = 1;
    outit = 0;
    resvec(it) = norm(r)/bnorm;
    x = x0;

    while(resvec(it) > tol) && (outit < maxit)
        outit = outit+1;
        rnorm = norm(r);
        [x,r,p,resvec_inner] = is_iter_gmres_svie(A,x,r,restart,L1,pL_L,L2,Mbprec,N_sie,tol*bnorm,rnorm);
        resvec(it+1:it+p) = resvec_inner/bnorm;
        it = it + p;
    end

    resvec = gather(resvec(1:it));
    iter = gather([p outit]);
    x = gather(x);
    clear r L1 L2 Mbprec;

end
