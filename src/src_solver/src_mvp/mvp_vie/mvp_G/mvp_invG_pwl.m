function T = mvp_invG_pwl(T,res)

    G = [1/res^3;1/(res^3/12);1/(res^3/12);1/(res^3/12);1/res^3;1/(res^3/12);1/(res^3/12);1/(res^3/12);1/res^3;1/(res^3/12);1/(res^3/12);1/(res^3/12)];

    for i = 1:12
        T(:,:,:,i) = G(i).*T(:,:,:,i);
    end
    
end