function [y] = mvp_coupling_Tucker_pwc(A,x,n1,n2,n3,m,r)
    
    y1 = zeros(n2*n1,n3,'like',x);
    y2 = zeros(n2*n1,n3,'like',x);
    y3 = zeros(n2*n1,n3,'like',x);
    
    for i = 1:m

        % Component 1
        t1 = A(i,1).G*x(i);        
        t2 = A(i,1).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,1).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,1).U3;
        y1  = y1 + t8;

        % Component 2
        t1 = A(i,2).G*x(i);        
        t2 = A(i,2).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,2).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,2).U3;
        y2 = y2 + t8;

        % Component 3
        t1 = A(i,3).G*x(i);        
        t2 = A(i,3).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,3).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,3).U3;
        y3 = y3 + t8;

    end
    
    y1 = reshape(y1,n2,n1,n3);
    y1 = permute(y1,[2,1,3]);
    y1 = y1(:);
    y2 = reshape(y2,n2,n1,n3);
    y2 = permute(y2,[2,1,3]);
    y2 = y2(:);
    y3 = reshape(y3,n2,n1,n3);
    y3 = permute(y3,[2,1,3]);
    y3 = y3(:);
    y = [y1;y2;y3];
    
end