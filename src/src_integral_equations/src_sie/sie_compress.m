function [Z_far_near] = sie_compress(Z,tol,comp_method)

    if strcmp(comp_method,'aca')
        [Z_far_near.U,Z_far_near.V,~,~] = aca(Z,tol);
    elseif strcmp(comp_method,'svd')
        [U,S,V] = svd(Z);
        [Z_far_near.U,S,V] = compress_SVD(U,S,V,tol);
        V = S*V';
        Z_far_near.V = V';
    end
   
end