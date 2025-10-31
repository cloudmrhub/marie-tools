function y = mvp_coupling_cross_TT_pwc(A,x)
    
    y1 = mvp_TT_core(A(1).TT,x);
    y2 = mvp_TT_core(A(2).TT,x);
    y3 = mvp_TT_core(A(3).TT,x);
    
    y = gather([y1;y2;y3]);
    
end
