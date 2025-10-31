function[G,U1,U2,U3] = hosvd(A,tol)
            
    [n1,n2,n3] = size(A);

    A1 = reshape(A,n1,[]);
    A2 = reshape(permute(A,[2,1,3]),n2,[]);
    A3 = reshape(permute(A,[3,1,2]),n3,[]);

    [U1,S1] = svd(A1,'econ');
    [U2,S2] = svd(A2,'econ');
    [U3,S3] = svd(A3,'econ');
    
    if nargin == 2
        r1 = find(diag(S1)<=tol*S1(1,1)/sqrt(3) & diag(S1)~=0,1);
        if isempty(r1)
            r1 = n1;
        end
        r2 = find(diag(S2)<=tol*S2(1,1)/sqrt(3) & diag(S2)~=0,1);
        if isempty(r2)
            r2 = n2;
        end
        r3 = find(diag(S3)<=tol*S3(1,1)/sqrt(3) & diag(S3)~=0,1);
        if isempty(r3)
            r3 = n3;
        end
    else
        r1 = rank(A1);
        r2 = rank(A2);
        r3 = rank(A3);
    end
    U1 = U1(:,1:r1);
    U2 = U2(:,1:r2);
    U3 = U3(:,1:r3);

    G = nmp(nmp(nmp(A,conj(U1),1),conj(U2),2),conj(U3),3);

end