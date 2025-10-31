function Jout = mvp_K_pwc_tucker(Jin,T)
    
    m1 = size(T(1,1).U1,2);
    m2 = size(T(1,1).U2,2);
    m3 = size(T(1,1).U3,2);
    
    [n1,n2,n3,pwx] = size(Jin);
    fJ  = zeros(m1,m2,m3,pwx,'like',Jin);
    KfJ = zeros(m1,m2,m3,'like',Jin);
    Jout  = zeros(n1,n2,n3,pwx,'like',Jin);

    pq = [0 -3 +2;+3 0 -1;-2 +1 0];

    for i=1:pwx
        fJ(:,:,:,i) = fftn(Jin(:,:,:,i),[m1,m2,m3]);
    end

    q_v = 1:3;

    for p = 1:3

        for q = q_v(1:end ~=p)

            j = abs(pq(p,q));

            KfJ = KfJ + sign(pq(p,q)) * nmp(nmp(nmp(T(j).G,T(j).U1,1),T(j).U2,2),T(j).U3,3) .* squeeze(fJ(:,:,:,q));

        end

        KfJ = ifftn(KfJ);
        Jout(:,:,:,p) = KfJ(1:n1,1:n2,1:n3);
        KfJ = zeros(m1,m2,m3,'like',Jin);
            
    end    

end