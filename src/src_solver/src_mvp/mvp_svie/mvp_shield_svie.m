function Jout = mvp_shield_svie(Jcb, ZscU, ZscV, Zbs, Zbc, Zss, Zcc, Mcr_inv, func, fN, emc, max_rank_shield_sie, max_rank_sie, N_shield_sie, N_sie, res, n1, n2, n3, ql, idx_solution)
                                     
    %% Currents
    Js               = Jcb(1:N_shield_sie,1);
    Jc               = Jcb(N_shield_sie+1:N_shield_sie+N_sie,1);
    Jb               = zeros(n1,n2,n3,ql,'like',Jcb);
    Jb(idx_solution) = Jcb(N_sie+1:end,1);
    
    %% VSIE solve
    if ( isa(Zss,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zss,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_ss = Zss*Js;
    else
        Jout_ss = gpuArray(Zss*gather(Js));
    end
    Jout_sc = ZscU*(ZscV*Jc);
    Jout_sb = gpuArray(func.mvp_transp_Zbc(Zbs,gather(Jb(:)),n1,n2,n3,N_wie,max_rank_shield_sie));

    if ( isa(Zcc,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zcc,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_cc = Zcc*Jc;
    else
        Jout_cc = gpuArray(Zcc*gather(Jc));
    end
    Jout_cs = ZscV.'*(ZscU.'*Js);
    Jout_cb = gpuArray(func.mvp_transp_Zbc(Zbc,gather(Jb(:)),n1,n2,n3,N_sie,max_rank_sie));
    
    Jout_bc = -gpuArray(func.mvp_Zbc(Zbc,gather(Jc),n1,n2,n3,N_sie,max_rank_sie));
    Jout_bc = Jout_bc(idx_solution);
    Jout_bb = Mcr_inv.*func.mvp_G(Jb,res)-func.mvp_N(Jb,fN);
    Jout_bb = 1 / emc.ce * Jout_bb(idx_solution);
    
    %% Concatenate and Return
    Jout   = [Jout_ss + Jout_sc + Jout_sb;
              Jout_cs + Jout_cc + Jout_cb;
              Jout_bs + Jout_bc + Jout_bb];

end