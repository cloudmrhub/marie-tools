function f = kernels_Nop(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,l,lp,np,ker_type)

%% inputs and outputs 

% X,Y,Z,Xp,Yp,Zp : N^4 integration points for the cases of constant planes {X,Y,Z} U {X',Y',Z'} transformed for integration over a unit square-square.
%                  IMPORTANT NOTE: the input argument for the testing and basis function is now r-r_m and r'-r'_n respectively.
% k0             : wavenumber
% dx             : resolution
% r_m            : (3x1) centre of the observation voxel
% r_n            : (3x1) centre of the source voxel
% np             : (3x1) normal vector on the surface of the source voxel
% ker_type       : \in {1,2,3,4} indicating the reduced surface-surface kernel to be calculated

% f              : complex double value of one surface-surface kernel


%% calculate the integral value for one surface-surface kernel

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
        
        R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
        G = 1/4/pi * exp(-1i*k0*R) ./ R;

        f = T .* B .* G;
        
    case 2
        %% surface-surface kernel 2
        
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

        q2 = np(1)*ff(:,1) + np(2)*ff(:,2) + np(3)*ff(:,3);

        f = B .* q2;
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

        %  n' . F
        ff(:,1) = (X-Xp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
        ff(:,2) = (Y-Yp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
        ff(:,3) = (Z-Zp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);

        q2 = np(1)*ff(:,1) + np(2)*ff(:,2) + np(3)*ff(:,3);

        f = T .* q2;
        
    case 4
        %% surface-surface kernel 4
        
        % Green
        R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
        G = 1/4/pi * exp(-1i*k0*R) ./ R;
        G0 = 1/4/pi ./ R;
        
        f = (G - G0) / (1i*k0)^2;
        
end


end