function [y] = mvp_coupling_Tucker_transpose_pwc(A,x,n1,n2,n3,m,r)

    x = reshape(x,n1,n2,n3,3);
    x = permute(x,[2,1,3,4]);
    x1 = x(:,:,:,1);
    x2 = x(:,:,:,2);
    x3 = x(:,:,:,3);
    x1 = x1(:).';
    x2 = x2(:).';
    x3 = x3(:).';
    
    y = zeros(m,1,'like',x);
    for i = 1:m
            
        % Component 1
        t2 = A(i,1).U1*A(i,1).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,1).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_1 = t7*A(i,1).U3;
        
        % Component 2
        t2 = A(i,2).U1*A(i,2).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,2).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_2 = t7*A(i,2).U3;
        
        % Component 3
        t2 = A(i,3).U1*A(i,3).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,3).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_3 = t7*A(i,3).U3;

        y(i) = x1*t8_1(:)+x2*t8_2(:)+x3*t8_3(:);

    end
    y = y(:);
    
end