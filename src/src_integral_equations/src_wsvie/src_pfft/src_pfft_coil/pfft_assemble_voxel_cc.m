function [pfft_Z_cc] = pfft_assemble_voxel_cc(MREDM,pfft_PS,Zcc,N_coil) 

    dims           = MREDM.dimensions;   
    ce             = MREDM.emc.ce;
    N_vie          = size(dims.idxS,1)*dims.ql;
    N_pfft         = dims.pfft_n1*dims.pfft_n2*dims.pfft_n3;
    pfft_kernel_n1 = dims.pfft_kernel_n1;
    pfft_kernel_n2 = dims.pfft_kernel_n2;
    pfft_kernel_n3 = dims.pfft_kernel_n3;
    res            = dims.res;
    pFFT_kernel    = dims.pFFT_kernel;
    ql             = dims.ql; 
    N_near_3D      = pfft_kernel_n1^3;
    N_pfft_dofs    = N_pfft*ql;
    N_near_dofs    = N_near_3D*ql;
    Jcb            = [speye(N_coil); sparse(N_vie, N_coil)];
    exp_cell_near  = MREDM.pfft.cell_id;
    c2c_near_list  = MREDM.pfft.C2C_near_list;
    near_vox_list  = MREDM.pfft.near_vox_list;
    exp_center     = floor(pfft_kernel_n1 / 2) + mod(pFFT_kernel, 2);
    id1            = pfft_proj_find_expansion_cell(exp_center, exp_center, exp_center, pFFT_kernel, pfft_kernel_n1, pfft_kernel_n1, pfft_kernel_n1); 
    s0             = size(id1,2);
    idx_exp_near   = zeros(s0*ql,1);
    for i = 1:ql
        idx_exp_near((i-1)*s0+1 : i*s0,1) = id1.' + (i-1)*N_near_3D;
    end
    pfft_n = to_GPU(MREDM.operators.pfft_n,1);

    Z_voxel  = cell(c2c_near_list.Count,1);
    I_idx    = cell(c2c_near_list.Count,1);
    J_idx    = cell(c2c_near_list.Count,1);

    for i_coil = 1:c2c_near_list.Count

        exp_voxels  = exp_cell_near(:,i_coil);
        near_list   = near_vox_list(:,i_coil);
        idx_near_1D = find(near_list);
        c2c_near    = c2c_near_list(num2str(i_coil));
        
        s1 = size(exp_voxels,1);
        s2 = size(nonzeros(near_list),1);
        s3 = size(idx_near_1D,1);
        idx_exp_pfft   = zeros(s1*ql,1);
        near_list_dofs = zeros(s2*ql,1);
        idx_near_3D    = zeros(s3*ql,1);
        for i = 1:ql
            idx_exp_pfft((i-1)*s1+1 : i*s1,1)   = exp_voxels + (i-1)*N_pfft;
            near_list_dofs((i-1)*s2+1 : i*s2,1) = nonzeros(near_list) + (i-1)*N_pfft;
            idx_near_3D((i-1)*s3+1 : i*s3,1)    = idx_near_1D + (i-1)*N_near_3D;
        end

        Jc_pfft               = pfft_PS * Jcb(:,i_coil);
        Jc_near               = zeros(N_near_dofs,1);
        Jc_near(idx_exp_near) = Jc_pfft(idx_exp_pfft);

        Jc_near_gpu    = reshape(gpuArray(Jc_near),[pfft_kernel_n1 pfft_kernel_n2 pfft_kernel_n3 ql]);
        Jout_near_gpu  = MREDM.functions.mvp_N(Jc_near_gpu,pfft_n);
        G_J            = MREDM.functions.mvp_G(Jc_near_gpu,res);
        Ec_near        = gather(Jout_near_gpu - G_J);
        
        Ec_pfft                 = zeros(N_pfft_dofs,1);
        Ec_pfft(near_list_dofs) = Ec_near(idx_near_3D);
        Ecb                     = pfft_PS.' * Ec_pfft;

        Z_voxel{i_coil} = 1/ce *Ecb(c2c_near);
        I_idx{i_coil}   = c2c_near';
        J_idx{i_coil}   = double(i_coil)*ones(length(c2c_near),1);  
    end
    Z_voxel    = vertcat(Z_voxel{:});
    I_idx      = vertcat(I_idx{:});
    J_idx      = vertcat(J_idx{:});
    Z_voxel_cc = sparse(I_idx, J_idx, Z_voxel, N_coil, N_coil);

    Z_direct = cell(c2c_near_list.Count,1);
    I_idx    = cell(c2c_near_list.Count,1);
    J_idx    = cell(c2c_near_list.Count,1);

    for i_coil = 1:N_coil
        c2c_near         = c2c_near_list(num2str(i_coil));
        Z_direct{i_coil} = Zcc(i_coil, c2c_near).';
        I_idx{i_coil}    = double(i_coil)*ones(length(c2c_near),1); 
        J_idx{i_coil}    = c2c_near';  
    end
    Z_direct    = vertcat(Z_direct{:});
    I_idx       = vertcat(I_idx{:});
    J_idx       = vertcat(J_idx{:});
    Z_direct_cc = sparse(I_idx, J_idx, Z_direct, N_coil, N_coil);

    pfft_Z_cc = Z_direct_cc - Z_voxel_cc;

end