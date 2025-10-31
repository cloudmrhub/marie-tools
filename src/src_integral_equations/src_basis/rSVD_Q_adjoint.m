function [Q] = rSVD_Q_adjoint(fK,func,dims,m,n,tol,blocksize,id1,id2)

    nb1                      = dims.nb1;
    nb2                      = dims.nb2;
    nb3                      = dims.nb3;
    res                      = dims.res;
    ql                       = dims.ql;
    
    flagdone = 0;
    ct = 0;
    
    AGn = zeros(m,blocksize);
    
    while (flagdone == 0)
        
        ct  = ct+1;
        An  = randn(n,ct*blocksize) + 1i*randn(n,ct*blocksize);
        
        for ii = blocksize*(ct-1)+1:blocksize*ct
            x         = zeros(nb1,nb2,nb3,ql);
            x(id1)    = An(:,ii)/norm(An(:,ii));
            x         = func.mvp_invG(x,res);
            x         = gather(-func.mvp_herm_adj_K(gpuArray(x),fK));
            AGn(:,ii) = x(id2);
        end
        
        [Q,S,~] = svd(AGn, 'econ'); 
        
        flagdone = 0;
        
        r1 = find(diag(S)<=tol*S(1,1) & diag(S)~=0,1);
        if ~isnan(r1)
            flagdone = 1;
            Q = Q(:,1:r1);
        else
            AGn = [AGn zeros(m,blocksize)];
        end
       
    end
end
    
