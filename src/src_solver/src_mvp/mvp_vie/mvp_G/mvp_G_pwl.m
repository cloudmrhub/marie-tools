function T = mvp_G_pwl(T,res)

    [~,~,~,pwx] = size(T);

    G = [res^3;res^3/12;res^3/12;res^3/12;res^3;res^3/12;res^3/12;res^3/12;res^3;res^3/12;res^3/12;res^3/12];
    
    for i = 1:pwx
        T(:,:,:,i) = G(i).*T(:,:,:,i);
    end
    
end