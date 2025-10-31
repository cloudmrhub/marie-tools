function [y] = mvp_coupling_Tucker_pwl(A,x,n1,n2,n3,m,r)
    
    y1 = zeros(n2*n1,n3,'like',x);
    y2 = zeros(n2*n1,n3,'like',x);
    y3 = zeros(n2*n1,n3,'like',x);
    y4 = zeros(n2*n1,n3,'like',x);
    y5 = zeros(n2*n1,n3,'like',x);
    y6 = zeros(n2*n1,n3,'like',x);
    y7 = zeros(n2*n1,n3,'like',x);
    y8 = zeros(n2*n1,n3,'like',x);
    y9 = zeros(n2*n1,n3,'like',x);
    y10 = zeros(n2*n1,n3,'like',x);
    y11 = zeros(n2*n1,n3,'like',x);
    y12 = zeros(n2*n1,n3,'like',x);
    
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

        % Component 4
        t1 = A(i,4).G*x(i);        
        t2 = A(i,4).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,4).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,4).U3;
        y4 = y4 + t8;

        % Component 5
        t1 = A(i,5).G*x(i);        
        t2 = A(i,5).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,5).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,5).U3;
        y5 = y5 + t8;

        % Component 6
        t1 = A(i,6).G*x(i);        
        t2 = A(i,6).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,6).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,6).U3;
        y6 = y6 + t8;

        % Component 7
        t1 = A(i,7).G*x(i);        
        t2 = A(i,7).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,7).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,7).U3;
        y7 = y7 + t8;

        % Component 8
        t1 = A(i,8).G*x(i);        
        t2 = A(i,8).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,8).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,8).U3;
        y8 = y8 + t8;

        % Component 9
        t1 = A(i,9).G*x(i);        
        t2 = A(i,9).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,9).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,9).U3;
        y9 = y9 + t8;

        % Component 10
        t1 = A(i,10).G*x(i);        
        t2 = A(i,10).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,10).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,10).U3;
        y10 = y10 + t8;

        % Component 11
        t1 = A(i,11).G*x(i);        
        t2 = A(i,11).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,11).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,11).U3;
        y11 = y11 + t8;

        % Component 12
        t1 = A(i,12).G*x(i);        
        t2 = A(i,12).U1*t1; 
        t3 = reshape(t2,n1,r,r);
        t4 = permute(t3,[2,1,3]);
        t5 = reshape(t4,r,n1*r);
        t6 = A(i,12).U2*t5;
        t7 = reshape(t6,n2*n1,r);
        t8 = t7*A(i,12).U3;
        y12 = y12 + t8;

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
    y4 = reshape(y4,n2,n1,n3);
    y4 = permute(y4,[2,1,3]);
    y4 = y4(:);
    y5 = reshape(y5,n2,n1,n3);
    y5 = permute(y5,[2,1,3]);
    y5 = y5(:);
    y6 = reshape(y6,n2,n1,n3);
    y6 = permute(y6,[2,1,3]);
    y6 = y6(:);
    y7 = reshape(y7,n2,n1,n3);
    y7 = permute(y7,[2,1,3]);
    y7 = y7(:);
    y8 = reshape(y8,n2,n1,n3);
    y8 = permute(y8,[2,1,3]);
    y8 = y8(:);
    y9 = reshape(y9,n2,n1,n3);
    y9 = permute(y9,[2,1,3]);
    y9 = y9(:);
    y10 = reshape(y10,n2,n1,n3);
    y10 = permute(y10,[2,1,3]);
    y10 = y10(:);
    y11 = reshape(y11,n2,n1,n3);
    y11 = permute(y11,[2,1,3]);
    y11 = y11(:);
    y12 = reshape(y12,n2,n1,n3);
    y12 = permute(y12,[2,1,3]);
    y12 = y12(:);
    y = [y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12];
    
end