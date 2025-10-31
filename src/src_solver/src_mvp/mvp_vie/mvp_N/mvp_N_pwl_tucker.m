function Jout = mvp_N_pwl_tucker(Jin,T)

    m1 = size(T(1,1).U1,2);
    m2 = size(T(1,1).U2,2);
    m3 = size(T(1,1).U3,2);
    
    [n1,n2,n3,pwx] = size(Jin);
    fJ  = zeros(m1,m2,m3,pwx,'like',Jin);
    NfJ = zeros(m1,m2,m3,'like',Jin);
    Jout  = zeros(n1,n2,n3,pwx,'like',Jin);

    pq = [1 2 3;2 4 5;3 5 6];
    llp = [1 5 6 7;-5 2 8 9;-6 8 3 10;-7 9 10 4];

    for i=1:pwx
        fJ(:,:,:,i) = fftn(Jin(:,:,:,i),[m1,m2,m3]);
    end

    for p = 1:3
        for l = 1:4

            left = 4*(p-1) + l;

            for q = 1:3
                for lp = 1:4

                    i = abs(llp(l,lp));
                    j = pq(p,q);

                    right = 4*(q-1) + lp;
                    NfJ  =  NfJ + sign(llp(l,lp)) * nmp(nmp(nmp(T(i,j).G,T(i,j).U1,1),T(i,j).U2,2),T(i,j).U3,3) .* squeeze(fJ(:,:,:,right));

                end
            end

            NfJ = ifftn(NfJ);
            Jout(:,:,:,left) = NfJ(1:n1,1:n2,1:n3);
            NfJ = zeros(m1,m2,m3,'like',Jin);

        end
    end
        
end
