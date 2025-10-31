function Vout = mvp_svie_pfft_tt(Jcb, pfft_PS, pfft_Z_cc, pfft_Z_bc_N, tt_Zbs_N, Zss, ZswcU, ZswcV, ext_Mc_inv, func, fN, Multiplier_1, Multiplier_2, emc, N_shield, N_coil, res, idx_solution, pfft_dims, nv, ql)
    
    %% Currents
    Jb                  = reshape(pfft_PS * [zeros(N_coil,1,'like',Jcb); Jcb(N_shield+N_coil+1:end,1)], pfft_dims);
    Jtot                = reshape(pfft_PS * Jcb(N_shield+1:end,1), pfft_dims);
    Jb_tt               = zeros(nv*ql,1,'like',Jcb);
    Jb_tt(idx_solution) = Jcb(N_shield+N_coil+1:end);
    
    %% VSIE solve
    Jout_pfft1   = func.mvp_N(Jtot,fN) - func.mvp_G(Jtot,res);
    Jout_pfft1   = 1/emc.ce * Multiplier_1.* (pfft_PS.' * Jout_pfft1(:));
    Jout_pfft2   = func.mvp_G(ext_Mc_inv.*Jb,res);
    Jout_pfft2   = 1/emc.ce * Multiplier_2.* (pfft_PS.' * Jout_pfft2(:));
    
    Jout_pfff3 = [pfft_Z_cc * Jcb(N_shield+1:N_shield+N_coil) + pfft_Z_bc_N.' * Jcb(N_shield+N_coil+1:end); -pfft_Z_bc_N * Jcb(N_shield+1:N_shield+N_coil)];
    
    Jout_tt_pfft = -func.mvp_Zbs(tt_Zbs_N,Jcb(1:N_shield));
    Jout_tt_pfft = [ZswcV.'*(ZswcU.'*Jcb(1:N_shield));Jout_tt_pfft(idx_solution)];
    
    if ( isa(Zss,'gpuArray') && isa(Jcb,'gpuArray') ) || ( ~isa(Zss,'gpuArray') && ~isa(Jcb,'gpuArray') )
        Jout_tt1 = Zss*Jcb(1:N_shield);
    else
        Jout_tt1 = gpuArray(Zss*gather(Jcb(1:N_shield)));
    end
    Jout_tt2     = ZswcU*(ZswcV*Jcb(N_shield+1:N_coil+N_shield)); 
    Jout_tt3     = func.mvp_transp_Zbs(tt_Zbs_N,Jb_tt,nv);
    
    %% Concatenate and Return 
    Vout = [Jout_tt1+Jout_tt2+Jout_tt3;Jout_pfft1+Jout_pfft2+Jout_pfff3+Jout_tt_pfft];   

end