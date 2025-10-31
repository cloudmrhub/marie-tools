function Jout = mvp_vie(Jb, EP, func, dim, fN)

    %% Dimensions
    res          = dim.res;
    n1           = dim.n1;
    n2           = dim.n2;
    n3           = dim.n3;
    ql           = dim.ql;
    idx_solution = dim.idx_solution;
    
    %% Currents
    Jp               = zeros(n1,n2,n3,ql,'like',Jb);
    Jp(idx_solution) = Jb;
    
    %% VIE solve
    N_J  = func.mvp_N(Jp,fN);
    Jout = Jp - (1./(EP.Mr).* EP.Mc).* func.mvp_invG(N_J,res);
    Jout = Jout(idx_solution);

end