function [Zbc_fft_N,Zbc_fft_K] = pfft_assemble_voxel_bc(MREDM,pfft_PS,N_coil)

    ce                 = MREDM.emc.ce;
    dims               = MREDM.dimensions;

    ql                 = dims.ql;
    pfft_idxS          = dims.pfft_idxS;
    pfft_kernel_n1     = dims.pfft_kernel_n1;
    pfft_kernel_n2     = dims.pfft_kernel_n2;
    pfft_kernel_n3     = dims.pfft_kernel_n3;
    pFFT_3Dkernel      = dims.pFFT_3Dkernel;
    N_pFFT_3Dkernel_ql = ql * pFFT_3Dkernel;
    N_near             = pfft_kernel_n1^3;
    N_near_ql          = ql * N_near;
    N_pfft             = dims.pfft_n1*dims.pfft_n2*dims.pfft_n3;
    N_b                = dims.N_b;
    N_scat             = dims.N_scat;
    res                = dims.res;
    cell_id            = MREDM.pfft.cell_id;
    S1_d               = MREDM.pfft.index_exp2near.S_1d;
    pwx_near_list      = MREDM.pfft.near_vox_list;
    C2B_near_list      = MREDM.pfft.C2B_near_list;
    keys               = C2B_near_list.keys;    
    idx_exp            = repmat(S1_d, ql,1) + kron(N_near * (0:ql-1).', ones(size(S1_d)));
    exp2pfft_dofs_idx  = repmat(cell_id, ql, 1) + kron(N_pfft * (0:ql-1).', ones(pFFT_3Dkernel, N_coil));
    
    Zbc_val_N        = cell(length(keys),1);
    Zbc_val_K        = cell(length(keys),1);
    I_idx            = Zbc_val_N;
    J_idx            = I_idx;
    
    Nnear_exp        = zeros(N_near_ql, N_pFFT_3Dkernel_ql);
    Knear_exp        = zeros(N_near_ql, N_pFFT_3Dkernel_ql);

    for i = 1:length(keys)
        i_sie = str2double(keys{i});
        exp_idx = exp2pfft_dofs_idx(:,i_sie);
        P_sie(:,i) = pfft_PS(exp_idx, i_sie);
    end
    
    try
        gpu_flag = 0;
        pfft_n = to_GPU(MREDM.operators.pfft_n,gpu_flag+1);
        pfft_k = to_GPU(MREDM.operators.pfft_k,gpu_flag+1);
        if ~isempty(keys)
            I                                 = eye(N_pFFT_3Dkernel_ql);
            Jin                               = zeros(N_near_ql, N_pFFT_3Dkernel_ql);
            Jin(idx_exp,1:N_pFFT_3Dkernel_ql) = I(:,1:N_pFFT_3Dkernel_ql);
            Jin                               = to_GPU(Jin,gpu_flag);
            for j = 1:N_pFFT_3Dkernel_ql
                Jin_near_gpu   = reshape(Jin(:,j),pfft_kernel_n1,pfft_kernel_n2,pfft_kernel_n3,ql);
                Jout_near_gpu  = gather(MREDM.functions.mvp_N(Jin_near_gpu,pfft_n));
                G_J            = gather(MREDM.functions.mvp_G(Jin_near_gpu,res));
                Nnear_exp(:,j) = reshape(Jout_near_gpu - G_J,[],1);
                Jout_near_gpu  = gather(MREDM.functions.mvp_K(Jin_near_gpu,pfft_k));
                Knear_exp(:,j) = reshape(Jout_near_gpu,[],1);
            end
            E_near = Nnear_exp * P_sie;
            H_near = Knear_exp * P_sie;
        end
    catch
        gpu_flag = 5;
        pfft_n = to_GPU(MREDM.operators.pfft_n,gpu_flag+1);
        pfft_k = to_GPU(MREDM.operators.pfft_k,gpu_flag+1);
        if ~isempty(keys)
            I                                 = eye(N_pFFT_3Dkernel_ql);
            Jin                               = zeros(N_near_ql, N_pFFT_3Dkernel_ql);
            Jin(idx_exp,1:N_pFFT_3Dkernel_ql) = I(:,1:N_pFFT_3Dkernel_ql);
            Jin                               = to_GPU(Jin,gpu_flag);
            for j = 1:N_pFFT_3Dkernel_ql
                Jin_near_gpu   = reshape(Jin(:,j),pfft_kernel_n1,pfft_kernel_n2,pfft_kernel_n3,ql);
                Jout_near_gpu  = gather(MREDM.functions.mvp_N(Jin_near_gpu,pfft_n));
                G_J            = gather(MREDM.functions.mvp_G(Jin_near_gpu,res));
                Nnear_exp(:,j) = reshape(Jout_near_gpu - G_J,[],1);
                Jout_near_gpu  = gather(MREDM.functions.mvp_K(Jin_near_gpu,pfft_k));
                Knear_exp(:,j) = reshape(Jout_near_gpu,[],1);
            end
            E_near = Nnear_exp * P_sie;
            H_near = Knear_exp * P_sie;
        end
    end

    for i = 1:length(keys)

        i_sie                = str2double(keys{i});
        near_scat_vox_pfft   = C2B_near_list(keys{i});
        near_vox_pfft        = pwx_near_list(:,i_sie);
        near_scat_vox_local  = find(ismember(near_vox_pfft, near_scat_vox_pfft)); 
        index_near_scat2dofs = find(ismember(pfft_idxS, near_scat_vox_pfft));
        index_near_scat_loc  = repmat(near_scat_vox_local, ql,1) + kron(N_near * (0:ql-1).', ones(size(near_scat_vox_local)));
        index_near_scat_dofs = repmat(index_near_scat2dofs, ql,1) + kron(N_scat * (0:ql-1).', ones(size(index_near_scat2dofs)));
        
        Zbc_val_N{i} = 1/ce*E_near(index_near_scat_loc,i);   
        Zbc_val_K{i} = H_near(index_near_scat_loc,i);      
        I_idx{i}     = index_near_scat_dofs;
        J_idx{i}     = i_sie * ones(size(index_near_scat_loc));

    end

    Zbc_val_N = vertcat(Zbc_val_N{:});
    Zbc_val_K = vertcat(Zbc_val_K{:});
    I_idx     = vertcat(I_idx{:});
    J_idx     = vertcat(J_idx{:});
    Zbc_fft_N = sparse(I_idx, J_idx, Zbc_val_N, N_b, N_coil);
    Zbc_fft_K = sparse(I_idx, J_idx, Zbc_val_K, N_b, N_coil);

end
