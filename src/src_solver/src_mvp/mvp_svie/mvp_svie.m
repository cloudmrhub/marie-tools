function Jout = mvp_svie(Jcb, Zbc, Zcc, Mcr_inv, func, fN, emc, max_rank_sie, N_sie, res, n1, n2, n3, ql, idx_solution)

    %% Currents
    Jc               = Jcb(1:N_sie,1);
    Jb               = zeros(n1,n2,n3,ql,'like',Jcb);
    Jb(idx_solution) = Jcb(N_sie+1:end,1);
    
    %% VSIE solve
    if ( isa(Zcc,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zcc,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_1 = Zcc*Jc;
    else
        Jout_1 = gpuArray(Zcc*gather(Jc));
    end
    Jout_2 = gpuArray(func.mvp_transp_Zbc(Zbc,gather(Jb(:)),n1,n2,n3,N_sie,max_rank_sie));
    Jout_3 = -gpuArray(func.mvp_Zbc(Zbc,gather(Jc),n1,n2,n3,N_sie,max_rank_sie));
    
    Jout_4 = Mcr_inv.*func.mvp_G(Jb,res)-func.mvp_N(Jb,fN);
    Jout_4 = 1 / emc.ce * Jout_4(idx_solution);
    
    %% Concatenate and Return
    Jout   = [Jout_1+Jout_2;Jout_3(idx_solution)+Jout_4];

end