function Vout = mvp_svie_pfft(Jcb, pfft_PS, pfft_Z_cc, pfft_Z_bc_N, ext_Mc_inv, func, fN, Multiplier_1, Multiplier_2, emc, N_sie, res, pfft_dims)
    
    %% Currents
    Jb   = reshape(pfft_PS * [zeros(N_sie,1,'like',Jcb); Jcb(N_sie+1:end,1)], pfft_dims);
    Jtot = reshape(pfft_PS * Jcb, pfft_dims);
    
    %% VSIE solve
    Jout_1 = func.mvp_N(Jtot,fN) - func.mvp_G(Jtot,res);
    Jout_1 = 1/emc.ce * Multiplier_1.* (pfft_PS.' * Jout_1(:));
    
    Jout_2 = func.mvp_G(ext_Mc_inv.*Jb,res);
    Jout_2 = 1/emc.ce * Multiplier_2.* (pfft_PS.' * Jout_2(:));
    
    Jout_3 = [pfft_Z_cc * Jcb(1:N_sie) + pfft_Z_bc_N.' * Jcb(N_sie+1:end); -pfft_Z_bc_N * Jcb(1:N_sie)];
    
    %% Concatenate and Return
    Vout = Jout_1+Jout_2+Jout_3;

end