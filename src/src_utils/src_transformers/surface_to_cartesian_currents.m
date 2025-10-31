function [Cart_p, j_cart] = surface_to_cartesian_currents(Jc,coil)

    coil_lumped_elements = geo_scoil_lumped_elements(coil,0);
    coil = geo_scoil(coil,coil_lumped_elements);
    
    first_node  = [3 1 2];
    second_node = [2 3 1];

    j_cart = zeros(size(coil.elem,2),3);
    for jj = 1:size(coil.elem,2)
        rp = coil.Ct(:,jj);
        r1 = coil.node(:,coil.elem(1,jj));
        r2 = coil.node(:,coil.elem(2,jj));
        r3 = coil.node(:,coil.elem(3,jj));
        Ae = Triangle_area(r1,r2,r3);
        Ae1 = Triangle_area(rp,r2,r3);
        Ae2 = Triangle_area(r1,rp,r3);
        Ae3 = Triangle_area(r1,r2,rp);
        zs = [Ae1 Ae2 Ae3]./Ae;
        ed = abs(coil.etod(:,jj));
        lo = [r2(1)-r3(1) r3(1)-r1(1) r1(1)-r2(1); r2(2)-r3(2) r3(2)-r1(2) r1(2)-r2(2); r2(3)-r3(3) r3(3)-r1(3) r1(3)-r2(3)];
        for id = 1:3
            if coil.index(ed(id))~=0
                i_1 = second_node(id);
                i_2 = first_node(id);
                L_i = sqrt(lo(1,id)^2 + lo(2,id)^2 + lo(3,id)^2);
                j_cart(jj,:) = j_cart(jj,:) + (L_i/(2*Ae)) * (zs(i_2)*lo(:,i_1)' - zs(i_1)*lo(:,i_2)') * sign(coil.etod(id,jj)) * Jc(coil.index(ed(id)),1);
            end
        end
    end
    Cart_p = [coil.Ct(1,:)' coil.Ct(2,:)' coil.Ct(3,:)'];
    
end