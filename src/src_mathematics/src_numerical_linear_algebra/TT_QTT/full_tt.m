function [a] = full_tt(tt)

    d=tt.d; 
    n=tt.n; 
    ps=tt.ps; 
    core=tt.core; 
    r=tt.r;
    
    a=core(ps(1):ps(2)-1);
    
    for i=2:d
      cr=core(ps(i):ps(i+1)-1);
      cr=reshape(cr,[r(i),n(i)*r(i+1)]);
      a=reshape(a,[numel(a)/r(i),r(i)]);
      a=a*cr;
    end

end
