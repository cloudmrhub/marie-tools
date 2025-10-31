function y = mvp_coupling_cross_TT_pwl(A,x)
    
    y1 = mvp_TT_core(A(1).TT,x);
    y2 = mvp_TT_core(A(2).TT,x);
    y3 = mvp_TT_core(A(3).TT,x);
    y4 = mvp_TT_core(A(4).TT,x);
    y5 = mvp_TT_core(A(5).TT,x);
    y6 = mvp_TT_core(A(6).TT,x);
    y7 = mvp_TT_core(A(7).TT,x);
    y8 = mvp_TT_core(A(8).TT,x);
    y9 = mvp_TT_core(A(9).TT,x);
    y10 = mvp_TT_core(A(10).TT,x);
    y11 = mvp_TT_core(A(11).TT,x);
    y12 = mvp_TT_core(A(12).TT,x);
    
    y = gather([y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12]);
    
end
