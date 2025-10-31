function [a] = mvp_TT_core(tt,x)

    n=tt.n; 
    ps=tt.ps; 
    r=tt.r;

    cr = tt.core(ps(4):ps(5)-1);
    cr = reshape(cr,[r(4),n(4)*r(5)]);
    a  = cr*x;
    
    cr = tt.core(ps(3):ps(4)-1);
    cr = reshape(cr,[r(3)*n(3),r(4)]);
    a  = cr*a;
    a  = reshape(a,[r(3),n(3)]);
    
    cr = tt.core(ps(2):ps(3)-1);
    cr = reshape(cr,[r(2)*n(2),r(3)]);
    a  = cr*a;
    a  = reshape(a,[r(2),n(2)*n(3)]);
    
    cr = tt.core(ps(1):ps(2)-1);
    cr = reshape(cr,[r(1)*n(1),r(2)]);
    a  = cr*a;
    a  = a(:);
    
end
