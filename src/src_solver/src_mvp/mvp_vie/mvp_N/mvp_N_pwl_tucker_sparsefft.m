function Jout = mvp_N_pwl_tucker_sparsefft(Jin,T)

    m1             = size(T(1,1).U1,2);
    m2             = size(T(1,1).U2,2);
    m3             = size(T(1,1).U3,2);
    [n1,n2,n3,pwx] = size(Jin);
    pq             = [1 2 3;2 4 5;3 5 6];
    llp            = [1 5 6 7;-5 2 8 9;-6 8 3 10;-7 9 10 4];
    
    fJ             = zeros(m1,m2,n3,pwx,'like',Jin);
    NfJ            = zeros(m1,m2,n3,'like',Jin);
    Jout           = zeros(n1,n2,n3,pwx,'like',Jin);
    
    for i=1:pwx
        fJ(:,:,:,i) = fft(fft(Jin(:,:,:,i),m1,1),m2,2);
    end

    b(1).a = 1:n1;
    b(1).b = 1:n2;
    b(2).a = 1:n1;
    b(2).b = n2+1:m2;
    b(3).a = n1+1:m1;
    b(3).b = 1:n2;
    b(4).a = n1+1:m1;
    b(4).b = n2+1:m2;

    for id = 1:4
    
        for p = 1:3
            for l = 1:4
                
                left = 4*(p-1) + l;
                
                for q = 1:3
                    for lp = 1:4
                            
                        i     = abs(llp(l,lp));
                        j     = pq(p,q);
                        right = 4*(q-1) + lp;
    
                        NfJ(b(id).a,b(id).b,:) = NfJ(b(id).a,b(id).b,:) + sign(llp(l,lp)) * nmp(nmp(nmp(T(i,j).G,T(i,j).U1(b(id).a,:),1),T(i,j).U2(b(id).b,:),2),T(i,j).U3,3) .* fft(fJ(b(id).a,b(id).b,:,right),m3,3);
                    
                    end
                end
                
                NfJ              = ifftn(NfJ);
                Jout(:,:,:,left) = NfJ(1:n1,1:n2,:);
                NfJ              = zeros(m1,m2,n3,'like',Jin);
                
            end
        end

    end

end