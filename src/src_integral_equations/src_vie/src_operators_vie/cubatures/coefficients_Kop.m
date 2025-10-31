function f = coefficients_Kop(dx,ko,n,np,p,q,l,lp,ker_type)

pp = zeros(3,1);
pp(p) = 1;

qq = zeros(3,1);
qq(q) = 1;


switch ker_type 
    
    case 1
        %% surface-surface kernel 1
        
        f = dot(cross(pp,qq),np);
        
        
    case 2
        %% surface-surface kernel 2
       
        switch l
            case 1
                h = [0 0 0];
            case 2
                h = [1/dx 0 0];
            case 3
                h = [0 1/dx 0];
            case 4
                h = [0 0 1/dx];
        end

        switch lp
            case 1
                hp = [0 0 0];
            case 2
                hp = [1/dx 0 0];
            case 3
                hp = [0 1/dx 0];
            case 4
                hp = [0 0 1/dx];
        end
        
        f = 1/(1i*ko)^2 * dot( np , dot(cross(qq,pp),hp)*h );

        
    case 3
        %% surface-surface kernel 3
       
        switch lp
            case 1
                hp = [0 0 0];
            case 2
                hp = [1/dx 0 0];
            case 3
                hp = [0 1/dx 0];
            case 4
                hp = [0 0 1/dx];
        end

       f = dot(n,np) * dot(hp,cross(pp,qq));

        
    case 4
        %% surface-surface kernel 4
        
        switch l
            case 1
                h = [0 0 0];
            case 2
                h = [1/dx 0 0];
            case 3
                h = [0 1/dx 0];
            case 4
                h = [0 0 1/dx];
        end


        a1 = cross(qq,[1 0 0]);
        a2 = cross(qq,[0 1 0]);
        a3 = cross(qq,[0 0 1]);

        b1 = h * dot(pp,[1 0 0]);
        b2 = h * dot(pp,[0 1 0]);
        b3 = h * dot(pp,[0 0 1]);

        f = dot(a1,np)*dot(b1,n) + dot(a2,np)*dot(b2,n) + dot(a3,np)*dot(b3,n);
        
end

end