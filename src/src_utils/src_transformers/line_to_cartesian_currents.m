function [j_cart,S_point] = line_to_cartesian_currents(Jw,filename)

    S = load(filename,'-ascii');
    loop_start       = find(squeeze(S(:,8))==1);
    loop_end         = find(squeeze(S(:,8))==2);
    F_point          = S(:,1:3);
    S_point          = S(:,4:6);
    T_point          = S(2:end,4:6);
    loop_start = loop_start(:).';  
    loop_end   = loop_end(:).';
    for i = 1:numel(loop_start)
        s = loop_start(i);
        e = loop_end(i);
        
        first_point = F_point(s:e,:);
        second_point = S_point(s:e,:);
        loop_length = length(F_point(s:e,:));

        c_s = s-loop_length*(i-1);
        c_e = e-loop_length*(i-1);

        F_point(s:e,:) = [first_point(c_e,:); first_point(c_s:c_e-1,:)];
        S_point(s:e,:) = first_point(c_s:c_e,:);
        T_point(s:e,:) = second_point(c_s:c_e,:);
    end
    
    Nep = size(T_point,1);
    Dl = vecnorm(S_point(1:Nep,:)-F_point(1:Nep,:),2,2);
    j_cart = zeros(Nep,3);

    for n = 1:Nep
        xs1 = F_point(n,1);
        ys1 = F_point(n,2);
        zs1 = F_point(n,3);
        xs2 = S_point(n,1); 
        ys2 = S_point(n,2);
        zs2 = S_point(n,3);
        t_left  = [(xs2-xs1), (ys2-ys1), (zs2-zs1)] / Dl(n);
        j_cart(n,:) = t_left * Jw(n);
    end
    j_cart = real(j_cart);
    
end