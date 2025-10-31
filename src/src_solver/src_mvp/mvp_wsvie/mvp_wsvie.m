function Jout = mvp_wsvie(Jcb, ZcwU, ZcwV, Zbw, Zbc, Zww, Zcc, Mcr_inv, func, fN, emc, max_rank_wie, max_rank_sie, N_wie, N_sie, res, n1, n2, n3, ql, idx_solution)

    %% Currents
    Jw               = Jcb(1:N_wie,1);
    Jc               = Jcb(1:N_wie+1:N_wie+N_sie,1);
    Jb               = zeros(n1,n2,n3,ql,'like',Jcb);
    Jb(idx_solution) = Jcb(N_wie+N_sie+1:end,1);
    
    %% VSWIE solve
    if ( isa(Zww,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zww,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_ww = Zww*Jw;
    else
        Jout_ww = gpuArray(Zww*gather(Jw));
    end
    Jout_wc = ZcwU*(ZcwV*Jc);
    Jout_wb = gpuArray(func.mvp_transp_Zbc(Zbw,gather(Jb(:)),n1,n2,n3,N_wie,max_rank_wie));

    if ( isa(Zcc,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zcc,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_cc = Zcc*Jc;
    else
        Jout_cc = gpuArray(Zcc*gather(Jc));
    end
    Jout_cw = ZcwV.'*(ZcwU.'*Jw);
    Jout_cb = gpuArray(func.mvp_transp_Zbc(Zbc,gather(Jb(:)),n1,n2,n3,N_sie,max_rank_sie));
    
    Jout_bb = Mcr_inv.*func.mvp_G(Jb,res)-func.mvp_N(Jb,fN);
    Jout_bb = 1 / emc.ce * Jout_bb(idx_solution);
    Jout_bw = -gpuArray(func.mvp_Zbc(Zbw,gather(Jw),n1,n2,n3,N_wie,max_rank_wie));
    Jout_bc = -gpuArray(func.mvp_Zbc(Zbc,gather(Jc),n1,n2,n3,N_sie,max_rank_sie));
    
    %% Concatenate and Return
    Jout = [Jout_ww+Jout_wc+Jout_wb;
            Jout_cc+Jout_cw++Jout_cb;
            Jout_bb+Jout_bw(idx_solution)+Jout_bc(idx_solution)];

end