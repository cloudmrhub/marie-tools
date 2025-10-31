function y = mvp_coupling_cross_TT_transpose_pwc(A,x,nv)
    
    y1 = mvp_transpose_TT_core(A(1,1).TT,x(1:nv,1));
    y2 = mvp_transpose_TT_core(A(1,2).TT,x(nv+1:2*nv,1));
    y3 = mvp_transpose_TT_core(A(1,3).TT,x(2*nv+1:3*nv,1));
    
    y = gather(y1.' + y2.' + y3.');
    
end
