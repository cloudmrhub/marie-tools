function Jout = mvp_svie_tt(Jcb, tt_Zbs_N, Zss, Mcr_inv, func, fN, emc, N_shield_sie, res, n1, n2, n3, ql, idx_solution)

    %% Currents
    Jc               = Jcb(1:N_shield_sie,1);
    Jb               = zeros(n1,n2,n3,ql,'like',Jcb);
    Jb(idx_solution) = Jcb(N_shield_sie+1:end,1);
    
    %% VSIE solve
    if ( isa(Zss,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zss,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_1 = Zss*Jc;
    else
        Jout_1 = gpuArray(Zss*gather(Jc));
    end
    Jout_2 = func.mvp_transp_Zbs(tt_Zbs_N,Jb(:),n1*n2*n3);
    Jout_3 = -func.mvp_Zbs(tt_Zbs_N,Jcb(1:N_shield_sie));
    
    Jout_4 = Mcr_inv.*func.mvp_G(Jb,res)-func.mvp_N(Jb,fN);
    Jout_4 = 1 / emc.ce * Jout_4(idx_solution);
    
    %% Concatenate and Return
    Jout   = [Jout_1+Jout_2;Jout_3(idx_solution)+Jout_4];

end