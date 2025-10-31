function [C] = nmp(A,B,mode)          
    switch mode
        case 1
            A = permute(A,[3,2,1]);
        case 2
            A = permute(A,[3,1,2]);
    end
    n1 = size(A,1);
    n2 = size(A,2);
    n3 = size(A,3);
    mn = size(B,2);
    A = reshape(A,n1*n2,n3);
    C = A*B;
    C = reshape(C,n1,n2,mn);
    switch mode
        case 1
            C = permute(C,[3,2,1]);
        case 2
            C = permute(C,[2,3,1]);    
    end
end

