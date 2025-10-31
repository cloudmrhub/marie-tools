function Jout = mvp_N_herm_adj_pwc_tucker(Jin,T)

    m1 = size(T(1,1).U1,2);
    m2 = size(T(1,1).U2,2);
    m3 = size(T(1,1).U3,2);
    
    [n1,n2,n3,pwx] = size(Jin);
    fJ  = zeros(m1,m2,m3,pwx,'like',Jin);
    NfJ = zeros(m1,m2,m3,'like',Jin);
    Jout  = zeros(n1,n2,n3,pwx,'like',Jin);

    pq = [1 2 3;2 4 5;3 5 6];

    for i=1:pwx
        fJ(:,:,:,i) = fftn(Jin(:,:,:,i),[m1,m2,m3]);
    end

    for p = 1:3

        for q = 1:3

            j = pq(p,q);
            NfJ  =  NfJ + conj(nmp(nmp(nmp(T(j).G,T(j).U1,1),T(j).U2,2),T(j).U3,3)) .* squeeze(fJ(:,:,:,q));

        end

        NfJ = ifftn(NfJ);
        Jout(:,:,:,p) = NfJ(1:n1,1:n2,1:n3);
        NfJ = zeros(m1,m2,m3,'like',Jin);

    end
        
end
