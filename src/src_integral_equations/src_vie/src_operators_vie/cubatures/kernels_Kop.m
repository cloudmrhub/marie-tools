function f = kernels_Kop(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,l,lp,n,ker_type)

switch ker_type
    case 1
        %% surface-surface kernel 1
        
        % testing function
        switch l
            case 1
                T =  1;
            case 2
                T = (X - r_m(1))/dx ;
            case 3
                T = (Y - r_m(2))/dx ;
            case 4
                T = (Z - r_m(3))/dx ;
        end
        
        % basis function
        switch lp
            case 1 
                B =  1;
            case 2
                B = (Xp - r_n(1))/dx ;
            case 3
                B = (Yp - r_n(2))/dx ;
            case 4
                B = (Zp - r_n(3))/dx ;
        end
        
        % Green
        R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
        G = 1/4/pi * exp(-1i*k0*R) ./ R;
        G0 = 1/4/pi ./ R;

        %  n' . F
        ff(:,1) = (X-Xp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
        ff(:,2) = (Y-Yp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
        ff(:,3) = (Z-Zp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);

        q2 = n(1)*ff(:,1) + n(2)*ff(:,2) + n(3)*ff(:,3);

        f = T .* B .* q2;
        
    case 2
        %% surface-surface kernel 2
        
        % Green
        R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
        G = 1/4/pi * exp(-1i*k0*R) ./ R;
        G0 = 1/4/pi ./ R;

        % F
        ff(:,1) = (X-Xp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
        ff(:,2) = (Y-Yp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
        ff(:,3) = (Z-Zp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);

        %F - RG0/2 
        gg(:,1) = (X-Xp)/2.*G0;
        gg(:,2) = (Y-Yp)/2.*G0;
        gg(:,3) = (Z-Zp)/2.*G0;
        hh = ff - gg;

        f = n(1)*hh(:,1) + n(2)*hh(:,2) + n(3)*hh(:,3);
        
    case 3
        %% surface-surface kernel 3
        
        % testing function
        switch l
            case 1
                T =  1;
            case 2
                T = (X - r_m(1))/dx ;
            case 3
                T = (Y - r_m(2))/dx ;
            case 4
                T = (Z - r_m(3))/dx ;
        end

        % Green
        R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
        G = 1/4/pi * exp(-1i*k0*R) ./ R;
        G0 = 1/4/pi ./ R;


        f = T .* (G-G0)/(1i*k0)^2;
        
    case 4
        %% surface-surface kernel 4
        
        % basis function
        switch lp
            case 1 
                B =  1;
            case 2
                B = (Xp - r_n(1))/dx ;
            case 3
                B = (Yp - r_n(2))/dx ;
            case 4
                B = (Zp - r_n(3))/dx ;
        end
        
        % Green
        R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
        G = 1/4/pi * exp(-1i*k0*R) ./ R;
        G0 = 1/4/pi ./ R;


        f = B .* (G-G0)/(1i*k0)^2;
        
end


end