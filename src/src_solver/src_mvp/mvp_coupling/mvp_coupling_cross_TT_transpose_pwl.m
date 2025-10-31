function y = mvp_coupling_cross_TT_transpose_pwl(A,x,nv)
    
    y1 = mvp_transpose_TT_core(A(1).TT,x(1:nv,1));
    y2 = mvp_transpose_TT_core(A(2).TT,x(nv+1:2*nv,1));
    y3 = mvp_transpose_TT_core(A(3).TT,x(2*nv+1:3*nv,1));
    y4 = mvp_transpose_TT_core(A(4).TT,x(3*nv+1:4*nv,1));
    
    y5 = mvp_transpose_TT_core(A(5).TT,x(4*nv+1:5*nv,1));
    y6 = mvp_transpose_TT_core(A(6).TT,x(5*nv+1:6*nv,1));
    y7 = mvp_transpose_TT_core(A(7).TT,x(6*nv+1:7*nv,1));
    y8 = mvp_transpose_TT_core(A(8).TT,x(7*nv+1:8*nv,1));
    
    y9 = mvp_transpose_TT_core(A(9).TT,x(8*nv+1:9*nv,1));
    y10 = mvp_transpose_TT_core(A(10).TT,x(9*nv+1:10*nv,1));
    y11 = mvp_transpose_TT_core(A(11).TT,x(10*nv+1:11*nv,1));
    y12 = mvp_transpose_TT_core(A(12).TT,x(11*nv+1:12*nv,1));
    
    y = gather(y1.' + y2.' + y3.' + y4.' + y5.' + y6.' + y7.' + y8.' + y9.' + y10.' + y11.' + y12.');
    
end
