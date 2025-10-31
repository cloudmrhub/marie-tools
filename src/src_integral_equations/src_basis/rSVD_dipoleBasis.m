function [BASIS] = rSVD_dipoleBasis(MREDM)

    tol_rSVD                 = MREDM.inputs.tol_rSVD;
    tol_DEIM                 = MREDM.inputs.tol_DEIM;
    blocksize                = MREDM.inputs.rSVD_blocksize;
    nb1                      = MREDM.dimensions.nb1;
    nb2                      = MREDM.dimensions.nb2;
    nb3                      = MREDM.dimensions.nb3;
    res                      = MREDM.dimensions.res;
    r                        = MREDM.dimensions.r;
    Nscat                    = MREDM.dimensions.N_scat;
    idxS                     = MREDM.dimensions.idxS;
    ql                       = MREDM.dimensions.ql;
    basis_idx_solution_solve = MREDM.dimensions.basis_idx_solution_solve;
    basis_idx_solution       = MREDM.dimensions.basis_idx_solution;
    m                        = length(MREDM.dimensions.basis_idx_solution_solve);
    n                        = length(MREDM.dimensions.basis_idx_solution);
    func                     = MREDM.functions;
    fK                       = to_GPU(MREDM.operators.basis_K,1);
    ce                       = MREDM.emc.ce;

    Q = rSVD_Q(fK,func,MREDM.dimensions,m,n,tol_DEIM,blocksize,basis_idx_solution,basis_idx_solution_solve);
    W = rSVD_Q_adjoint(fK,func,MREDM.dimensions,n,m,tol_DEIM,blocksize,basis_idx_solution_solve,basis_idx_solution);

    if size(Q,2) < size(W,2)
        Ub = zeros(n,size(Q,2));
        for i = 1:size(Q,2)
            x                           = zeros(nb1,nb2,nb3,ql);
            x(basis_idx_solution_solve) = Q(:,i);
            x                           = func.mvp_invG(x,res);
            x                           = gather(func.mvp_herm_adj_K(gpuArray(x),fK));
            Ub(:,i)                     = x(basis_idx_solution);
        end
        T = Ub'*W;
    else
        Ub = zeros(m,size(W,2));
        for i = 1:size(W,2)
            x                           = zeros(nb1,nb2,nb3,ql);
            x(basis_idx_solution)       = W(:,i);
            x                           = gather(func.mvp_K(gpuArray(x),fK));
            x                           = func.mvp_invG(x,res);
            Ub(:,i)                     = x(basis_idx_solution_solve);
        end
        T = Q'*Ub;
    end
    clear fK  % no longer needed

    [Ub, SK, VK] = svd(T,'econ');
    clear T % no longer needed

    r1 = find(diag(SK)<=tol_rSVD*SK(1,1) & diag(SK)~=0,1);

    BASIS.Ub = Q*Ub;
    clear Ub % no longer needed
    VK = W*VK;
    clear Q W % no longer needed

    [P,ndeim,xds,yds,zds] = rSVD_deim(BASIS.Ub,ql,Nscat,idxS,r);

    if ~isnan(r1)
        BASIS.Ub = BASIS.Ub(:,1:r1);
        SK = SK(1:r1,1:r1);
        VK = VK(:,1:r1);
    end

    fN = to_GPU(MREDM.operators.basis_N,1);
    BASIS.Ue = zeros(size(BASIS.Ub));
    for ii = 1:size(BASIS.Ue,2)
        x                           = zeros(nb1,nb2,nb3,ql);
        x(basis_idx_solution)       = VK(:,ii);
        x                           = 1/ce*gather(func.mvp_N(gpuArray(x),fN)) / SK(ii,ii);
        x                           = func.mvp_invG(x,res);
        BASIS.Ue(:,ii)              = x(basis_idx_solution_solve);
    end
    clear fN  % no longer needed

    X = (P.'*BASIS.Ue)\speye(ndeim,ndeim);

    BASIS.SK  = SK;
    BASIS.VK  = VK;
    BASIS.P   = P;
    BASIS.X   = X;
    BASIS.xds = xds;
    BASIS.yds = yds;
    BASIS.zds = zds;

    clear SK VK P X xds yds zds % no longer needed

    fprintf('\trank=%i, last singular value = %4.7f.\n',size(BASIS.Ub,2),BASIS.SK(end,end)/BASIS.SK(1,1));

end
