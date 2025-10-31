function[T] = assembly_fft_circ_tucker_pwl(A,MREDM,n1,n2,n3)

    [~,~,~,n4,n5] = size(A);
    tol = MREDM.inputs.tol_HOSVD;
    
    T = struct('G',repmat({[]}, n4, n5),'U1',repmat({[]}, n4, n5),'U2',repmat({[]}, n4, n5),'U3',repmat({[]}, n4, n5));

    if n5 == 6
        G_n1  = [+1 -1 -1 +1 +1 +1];
        G_n2  = [+1 -1 +1 +1 -1 +1];
        G_n3  = [+1 +1 -1 +1 -1 +1];
    elseif n5 == 3
        G_n1 = [-1 +1 +1];
        G_n2 = [+1 -1 +1];
        G_n3 = [+1 +1 -1];
    end
    bt_n1 = [+1 +1 +1 +1 -1 +1 +1 -1 -1 +1]; 
    bt_n2 = [+1 +1 +1 +1 +1 -1 +1 -1 +1 -1]; 
    bt_n3 = [+1 +1 +1 +1 +1 +1 -1 +1 -1 -1]; 
    
    for j=1:n4
        for k=1:n5
            
            [T(j,k).G,U1,U2,U3] = hosvd(A(:,:,:,j,k),tol);
            T(j,k).U1 = (fft([U1;zeros(1,size(U1,2));U1(n1:-1:2,:)*G_n1(k)*bt_n1(j)],[],1) ).';
            T(j,k).U2 = (fft([U2;zeros(1,size(U2,2));U2(n2:-1:2,:)*G_n2(k)*bt_n2(j)],[],1) ).';
            T(j,k).U3 = (fft([U3;zeros(1,size(U3,2));U3(n3:-1:2,:)*G_n3(k)*bt_n3(j)],[],1) ).';
            
        end
    end
    
end
