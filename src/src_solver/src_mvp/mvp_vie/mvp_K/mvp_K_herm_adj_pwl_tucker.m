function Jout = mvp_K_herm_adj_pwl_tucker(Jin,T)
    
    m1 = size(T(1,1).U1,2);
    m2 = size(T(1,1).U2,2);
    m3 = size(T(1,1).U3,2);
    
    [n1,n2,n3,pwx] = size(Jin);
    fJ  = zeros(m1,m2,m3,pwx,'like',Jin);
    KfJ = zeros(m1,m2,m3,'like',Jin);
    Jout  = zeros(n1,n2,n3,pwx,'like',Jin);

    pq = [0 -3 +2;+3 0 -1;-2 +1 0];
    llp = [1 5 6 7;-5 2 8 9;-6 8 3 10;-7 9 10 4];    
    slp = [1 -1 -1 -1;-1 1 1 1;-1 1 1 1 ;-1  1  1 1];

    for i=1:pwx
        fJ(:,:,:,i) = fftn(Jin(:,:,:,i),[m1,m2,m3]);
    end

    q_v = 1:3;

    for p = 1:3
        for l = 1:4
            
            left = 4*(p-1) + l;
            
            for q = q_v(1:end ~=p)
                for lp = 1:4
                        
                    i = abs(llp(l,lp));
                    j = abs(pq(p,q));
                
                    right = 4*(q-1) + lp;
                    KfJ = KfJ + sign(pq(p,q)) * sign(slp(l,lp)) * sign(llp(l,lp)) * conj(-nmp(nmp(nmp(T(i,j).G,T(i,j).U1,1),T(i,j).U2,2),T(i,j).U3,3)) .* fJ(:,:,:,right);
            
                end
            end
            
            KfJ = ifftn(KfJ);
            Jout(:,:,:,left) = KfJ(1:n1,1:n2,1:n3);
            KfJ = zeros(m1,m2,m3,'like',Jin);
            
        end
    end
    Jout = gather(Jout); 
end