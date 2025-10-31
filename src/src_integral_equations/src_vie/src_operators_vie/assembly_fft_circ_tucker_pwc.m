function[T] = assembly_fft_circ_tucker_pwc(A,MREDM,n1,n2,n3)

    [~,~,~,n4] = size(A);  
    tol = MREDM.inputs.tol_HOSVD;
    
    T = struct('G',repmat({[]}, n4,1),'U1',repmat({[]}, n4,1),'U2',repmat({[]}, n4,1),'U3',repmat({[]}, n4,1));

    if n4 == 6
        G_n1  = [+1 -1 -1 +1 +1 +1];
        G_n2  = [+1 -1 +1 +1 -1 +1];
        G_n3  = [+1 +1 -1 +1 -1 +1];
    elseif n4 == 3
        G_n1 = [-1 +1 +1];
        G_n2 = [+1 -1 +1];
        G_n3 = [+1 +1 -1];
    end
    
    for j=1:n4

        [T(j).G,U1,U2,U3] = hosvd(A(:,:,:,j),tol);

        T(j).U1 = (fft([U1;zeros(1,size(U1,2));U1(n1:-1:2,:)*G_n1(j)],[],1) ).';
        T(j).U2 = (fft([U2;zeros(1,size(U2,2));U2(n2:-1:2,:)*G_n2(j)],[],1) ).';
        T(j).U3 = (fft([U3;zeros(1,size(U3,2));U3(n3:-1:2,:)*G_n3(j)],[],1) ).';
            
    end
    
end
