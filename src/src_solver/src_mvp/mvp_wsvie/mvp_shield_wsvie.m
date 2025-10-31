function Jout = mvp_shield_wsvie(Jcb, ZswU, ZswV, ZscU, ZscV, ZwcU, ZwcV, Zbs, Zbw, Zbc, Zss, Zww, Zcc, Mcr_inv, func, fN, emc, max_rank_shield_sie, max_rank_wie, max_rank_sie, N_shield_sie, N_wie, N_sie, res, n1, n2, n3, ql, idx_solution)
    %% Currents
    Js               = Jcb(1:N_shield_sie,1);
    Jw               = Jcb(N_shield_sie+1:N_shield_sie+N_wie,1);
    Jc               = Jcb(N_shield_sie+N_wie+1:N_shield_sie+N_wie+N_sie,1);
    Jb               = zeros(n1,n2,n3,ql,'like',Jcb);
    Jb(idx_solution) = Jcb(N_wie+N_sie+1:end,1);
    
    %% VSWIE solve
    if ( isa(Zss,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zss,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_ss = Zss*Js;
    else
        Jout_ss = gpuArray(Zss*gather(Js));
    end
    Jout_sw = ZswU*(ZswV*Jw);
    Jout_sc = ZscU*(ZscV*Jc);
    Jout_sb = gpuArray(func.mvp_transp_Zbc(Zbs,gather(Jb(:)),n1,n2,n3,N_wie,max_rank_shield_sie));

    if ( isa(Zww,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zww,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_ww = Zww*Jw;
    else
        Jout_ww = gpuArray(Zww*gather(Jw));
    end
    Jout_wc = ZwcU*(ZwcV*Jc);
    Jout_ws = ZswV.'*(ZswV.'*Js);
    Jout_wb = gpuArray(func.mvp_transp_Zbc(Zbw,gather(Jb(:)),n1,n2,n3,N_wie,max_rank_wie));

    if ( isa(Zcc,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zcc,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_cc = Zcc*Jc;
    else
        Jout_cc = gpuArray(Zcc*gather(Jc));
    end
    Jout_cw = ZwcV.'*(ZwcU*Jw);
    Jout_cs = ZscV.'*(ZscU*Js);
    Jout_cb = gpuArray(func.mvp_transp_Zbc(Zbc,gather(Jb(:)),n1,n2,n3,N_sie,max_rank_sie));
    
    Jout_bb = Mcr_inv.*func.mvp_G(Jb,res)-func.mvp_N(Jb,fN);
    Jout_bb = 1 / emc.ce * Jout_bb(idx_solution);
    Jout_bs = -gpuArray(func.mvp_Zbc(Zbs,gather(Js),n1,n2,n3,N_wie,max_rank_shield_sie));
    Jout_bw = -gpuArray(func.mvp_Zbc(Zbw,gather(Jw),n1,n2,n3,N_wie,max_rank_wie));
    Jout_bw = Jout_bw(idx_solution);
    Jout_bc = -gpuArray(func.mvp_Zbc(Zbc,gather(Jc),n1,n2,n3,N_sie,max_rank_sie));
    Jout_bc = Jout_bc(idx_solution);
    
    %% Concatenate and Return
    Jout = [Jout_ss + Jout_sw + Jout_sc + Jout_sb;
            Jout_ws + Jout_ww + Jout_wc + Jout_wb;
            Jout_cs + Jout_cw + Jout_cc + Jout_cb;
            Jout_bs + Jout_bw + Jout_bc + Jout_bb];

end