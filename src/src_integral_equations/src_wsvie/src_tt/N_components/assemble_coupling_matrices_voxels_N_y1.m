function[Zbc_N] = assemble_coupling_matrices_voxels_N_y1(SIE_quads,VIE_quads,k0,n1,n2,n3,xd,yd,zd,res,rp,rn,r2,r3,ttnS)

    [ttnS,ii] = sortrows(ttnS,4);
    [ttn4,kafs] = unique(ttnS(:,4),'stable');
    
    kaf = zeros(size(kafs),'like',kafs);
    kbe = zeros(size(kafs),'like',kafs);
    
    kaf(end)     = size(ttnS,1);
    kaf(1:end-1) = kafs(2:end)-1;
    kbe(1)       = 1;
    kbe(2:end)   = kaf(1:end-1) + 1;
    
    Scoord = cell(length(ttn4),1);
    Zbc_Nt = cell(length(ttn4),1);
    
    for i = 1:length(ttn4)
        ttn1      = ttnS(kbe(i):kaf(i),1);
        ttn2      = ttnS(kbe(i):kaf(i),2);
        ttn3      = ttnS(kbe(i):kaf(i),3);
        ttN_v     = sub2ind([n1 n2 n3],ttn1,ttn2,ttn3);
        Scoord{i} = [xd(ttN_v) yd(ttN_v) zd(ttN_v)];
    end
    
    rpt = rp(:,ttn4);
    rnt = rn(:,ttn4);
    r2t = r2(:,ttn4);
    r3t = r3(:,ttn4);
       
    for i = 1:length(ttn4)
        Zbc_Nt{i} = Assemble_rwg_coupling_matrix_N_y1(Scoord{i}.',rpt(:,i),rnt(:,i),r2t(:,i),r3t(:,i),SIE_quads,VIE_quads,res,k0);
    end
    
    Zbc_N  = cell2mat(Zbc_Nt)*res^3;
    Zbc_N  = gpuArray(Zbc_N);
    Zbc_N(ii) = Zbc_N;
    
end