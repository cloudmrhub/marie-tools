function[Ae] = Triangle_area(r1,r2,r3)
    
    ls = [r2(1)-r3(1) r3(1)-r1(1) r1(1)-r2(1); r2(2)-r3(2) r3(2)-r1(2) r1(2)-r2(2); r2(3)-r3(3) r3(3)-r1(3) r1(3)-r2(3)];
    
    aa = norm(ls(:,1));
    bb = norm(ls(:,2));
    cc = norm(ls(:,3));
    
    s = (aa+bb+cc)/2;
    
    Ae = sqrt(s*(s-aa)*(s-bb)*(s-cc));
    
end