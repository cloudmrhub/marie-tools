function [y] = mvp_coupling_Tucker_transpose_pwl(A,x,n1,n2,n3,m,r)

    x = reshape(x,n1,n2,n3,12);
    x = permute(x,[2,1,3,4]);
    x1 = x(:,:,:,1);
    x2 = x(:,:,:,2);
    x3 = x(:,:,:,3);
    x4 = x(:,:,:,4);
    x5 = x(:,:,:,5);
    x6 = x(:,:,:,6);
    x7 = x(:,:,:,7);
    x8 = x(:,:,:,8);
    x9 = x(:,:,:,9);
    x10 = x(:,:,:,10);
    x11 = x(:,:,:,11);
    x12 = x(:,:,:,12);
    x1 = x1(:).';
    x2 = x2(:).';
    x3 = x3(:).';
    x4 = x4(:).';
    x5 = x5(:).';
    x6 = x6(:).';
    x7 = x7(:).';
    x8 = x8(:).';
    x9 = x9(:).';
    x10 = x10(:).';
    x11 = x11(:).';
    x12 = x12(:).';
    
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
        
        % Component 4
        t2 = A(i,4).U1*A(i,4).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,4).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_4 = t7*A(i,4).U3;
        
        % Component 5
        t2 = A(i,5).U1*A(i,5).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,5).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_5 = t7*A(i,5).U3;
        
        % Component 6
        t2 = A(i,6).U1*A(i,6).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,6).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_6 = t7*A(i,6).U3;
            
        % Component 7
        t2 = A(i,7).U1*A(i,7).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,7).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_7 = t7*A(i,7).U3;
        
        % Component 8
        t2 = A(i,8).U1*A(i,8).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,8).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_8 = t7*A(i,8).U3;
        
        % Component 9
        t2 = A(i,9).U1*A(i,9).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,9).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_9 = t7*A(i,9).U3;
        
        % Component 10
        t2 = A(i,10).U1*A(i,10).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,10).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_10 = t7*A(i,10).U3;
        
        % Component 11
        t2 = A(i,11).U1*A(i,11).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,11).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_11 = t7*A(i,11).U3;
        
        % Component 12
        t2 = A(i,12).U1*A(i,12).G; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,12).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8_12 = t7*A(i,12).U3;

        y(i) = x1*t8_1(:)+x2*t8_2(:)+x3*t8_3(:)+...
               x4*t8_4(:)+x5*t8_5(:)+x6*t8_6(:)+...
               x7*t8_7(:)+x8*t8_8(:)+x9*t8_9(:)+...
               x10*t8_10(:)+x11*t8_11(:)+x12*t8_12(:);

    end
    y = y(:);
    
end