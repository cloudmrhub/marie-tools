function [a] = mvp_transpose_TT_core(tt,x)

    n=tt.n; 
    ps=tt.ps; 
    r=tt.r;
    
    x = reshape(x,n(1),n(2),n(3));
    x = permute(x,[3 2 1]);
    x = reshape(x,n(3)*n(2),n(1));
    
    cr = tt.core(ps(1):ps(2)-1);
    cr = reshape(cr,[n(1),r(2)]);
    a  = x*cr;
    a  = reshape(a,[n(3),n(2),r(2)]);
    a  = permute(a,[1 3 2]);
    a  = reshape(a,n(3),r(2)*n(2));
    
    cr = tt.core(ps(2):ps(3)-1);
    cr = reshape(cr,[r(2)*n(2),r(3)]);
    a  = a*cr;
    a  = permute(a,[2 1]);
    a  = reshape(a,1,r(3)*n(3));
    
    cr = tt.core(ps(3):ps(4)-1);
    cr = reshape(cr,[r(3)*n(3),r(4)]);
    a  = a*cr;
    a  = reshape(a,1,r(4));
    
    cr = tt.core(ps(4):ps(5)-1);
    cr = reshape(cr,[r(4),n(4)]);
    a  = a*cr;

end