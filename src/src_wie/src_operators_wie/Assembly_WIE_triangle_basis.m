function [Z,Z_losses] = Assembly_WIE_triangle_basis(wire,a,order,emc)
    % first_p: first point of segment  (previous point) n-1
    % second_p: second point of segment (current point) n
    % third_p: third point of segment (next point) n+1
    % a: radius of the wire
    % emc: electromagnetic constants

    first_p  = wire.F_point;
    second_p = wire.S_point;
    third_p  = wire.T_point;

    % Wavenumber
    k0       = emc.k0;
    omega    = emc.omega;
    mu0      = emc.mu0;
    sc_const = 1i*omega*mu0/(4*pi);

    % Number of functions
    Nep = size(third_p,1);

    % Distances between points
    Dl = vecnorm(second_p(1:Nep,:)-first_p(1:Nep,:),2,2);
    Dr = vecnorm(third_p(1:Nep,:)-second_p(1:Nep,:),2,2);

    % MoM matrix
    Z = zeros(Nep,Nep);
    % Z_losses = zeros(Nep,Nep); 

    % Gauss quadrature integration
    [Weights,Xi] = gauss_1d(order);
    Weights = [Weights;Weights];

    % Loop around all points
    for n = 1:Nep

        % Source points
        xs1 = first_p(n,1);
        ys1 = first_p(n,2);
        zs1 = first_p(n,3);
        xs2 = second_p(n,1); 
        ys2 = second_p(n,2);
        zs2 = second_p(n,3);
        xs3 = third_p(n,1);
        ys3 = third_p(n,2);
        zs3 = third_p(n,3);

        % Define constant tangent vectors (the segments are straight)
        t_left  = [(xs2-xs1), (ys2-ys1), (zs2-zs1)] / Dl(n);
        t_right = [(xs3-xs2), (ys3-ys2), (zs3-zs2)] / Dr(n);
        tlp = zeros(3, 2*order);

        % Basis function
        fn  = zeros(2*order,1);
        dfn = zeros(2*order,1);

        % Source coordinates
        XS  = zeros(2*order,1); 
        YS  = zeros(2*order,1);  
        ZS  = zeros(2*order,1); 

        % Parametric integration Point
        lp  = zeros(2*order,1); 

        for p = 1:order

            % Left Segment
            lp(p)        = Dl(n)/2*(Xi(p) + 1);

            fn(p)        = lp(p)/Dl(n);
            dfn(p)       = 1/Dl(n);

            XS(p)        = xs1 + lp(p)/Dl(n)*(xs2-xs1);
            YS(p)        = ys1 + lp(p)/Dl(n)*(ys2-ys1);
            ZS(p)        = zs1 + lp(p)/Dl(n)*(zs2-zs1);
            tlp(:,p)     = t_left(:); % assign constant tangent vector for left segment


            % Right Segment
            lp(p+order)  = Dr(n)/2*(Xi(p) + 1);

            fn(p+order)  = (Dr(n)-lp(p+order))/Dr(n);
            dfn(p+order) = -1/Dr(n);

            XS(p+order)  = xs2 + lp(p+order)/Dr(n)*(xs3-xs2);
            YS(p+order)  = ys2 + lp(p+order)/Dr(n)*(ys3-ys2);
            ZS(p+order)  = zs2 + lp(p+order)/Dr(n)*(zs3-zs2);
            tlp(:,p+order) = t_right(:); % assign constant tangent vector for right segment

        end

        %% Self Term
        % Left-Left Segment
        Z_update_left_left = 0;
        for q = 1:order
            sqrt1              = sqrt(a^2 + (lp(q)-Dl(n))^2);
            sqrt2              = sqrt(a^2 + lp(q)^2);
            logterm            = log( abs(lp(q)+sqrt2) / abs(lp(q)-Dl(n)+sqrt1) );
            term               = 1i*k0*Dl(n)/2;
            S1                 = lp(q)/Dl(n)*logterm + 1/Dl(n)*sqrt1 - 1/Dl(n)*sqrt2 - term;
            S2                 = 1/Dl(n) * (logterm - 2*term);
            Z_update_left_left = Z_update_left_left + Weights(q)*(fn(q)*S1-1/(k0^2)*dfn(q)*S2);
        end
        Z_update_left_left = Dl(n)/2 * Z_update_left_left;

        % Left-Right Segment
        Z_update_left_right = 0;
        for p = 1:order
            for q = order+1:2*order
                T = fn(p)*fn(q)*dot(tlp(:,p),tlp(:,q)) -1/k0^2 * dfn(p)*dfn(q);
                R = sqrt((XS(p)-XS(q))^2+(YS(p)-YS(q))^2+(ZS(p)-ZS(q))^2 + a^2);
                G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                Z_update_left_right = Z_update_left_right + G*T;
            end
        end
        Z_update_left_right = Dl(n)/2 * Dr(n)/2 * Z_update_left_right;

        % Right-Left Segment interactions
        Z_update_right_left = 0;
        for p = order+1:2*order
            for q = 1:order
                T = fn(p)*fn(q)*dot(tlp(:,p),tlp(:,q)) -1/k0^2 * dfn(p)*dfn(q);
                R = sqrt((XS(p)-XS(q))^2+(YS(p)-YS(q))^2+(ZS(p)-ZS(q))^2 + a^2);
                G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                Z_update_right_left = Z_update_right_left + G*T;
            end
        end
        Z_update_right_left = Dr(n)/2 * Dl(n)/2 * Z_update_right_left;

        % Right-Right Segment
        Z_update_right_right = 0;
        for q = order+1:2*order
            sqrt1                = sqrt(a^2 + (lp(q)-Dr(n))^2 );
            sqrt2                = sqrt(a^2 + lp(q)^2 );
            logterm              = log( abs(lp(q)+sqrt2) / abs(lp(q)-Dr(n)+sqrt1) );
            term                 = 1i*k0*Dr(n)/2;
            S1                   = (1-lp(q)/Dr(n))*logterm - 1/Dr(n)*sqrt1 + 1/Dr(n)*sqrt2 - term;
            S2                   = -1/Dr(n) * (logterm - 2*term);
            Z_update_right_right = Z_update_right_right + Weights(q)*(fn(q)*S1-1/(k0^2)*dfn(q)*S2);
        end
        Z_update_right_right = Dr(n)/2 * Z_update_right_right;

        % Add all together
        Z(n,n) = Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right;

        %% Copper Losses Term
        % Left-Left Segment
        % Z_update_left_left = 0;
        % for p = 1:order
        %     for q = 1:order
        %         T = fn(p)*fn(q)*dot(tlp(:,p),tlp(:,q));
        %         Z_update_left_left = Z_update_left_left + Weights(p)*Weights(q)*T;
        %     end
        % end
        % % Z_update_left_left = (Dl(n)/2)^2 * Z_update_left_left;
        % 
        % % Left-Right Segment
        % Z_update_left_right = 0;
        % for p = 1:order
        %     for q = order+1:2*order
        %         T = fn(p)*fn(q)*dot(tlp(:,p),tlp(:,q));
        %         Z_update_left_right = Z_update_left_right + Weights(p)*Weights(q)*T;
        %     end
        % end
        % % Z_update_left_right = Dl(n)/2 * Dr(n)/2 * Z_update_left_right;
        % 
        % % Right-Left Segment interactions
        % Z_update_right_left = 0;
        % for p = order+1:2*order
        %     for q = 1:order
        %         T = fn(p)*fn(q)*dot(tlp(:,p),tlp(:,q));
        %         Z_update_right_left = Z_update_right_left + Weights(p)*Weights(q)*T;
        %     end
        % end
        % % Z_update_right_left = Dr(n)/2 * Dl(n)/2 * Z_update_right_left;
        % 
        % % Right-Right Segment
        % Z_update_right_right = 0;
        % for p = order+1:2*order
        %     for q = order+1:2*order
        %         T = fn(p)*fn(q)*dot(tlp(:,p),tlp(:,q));
        %         Z_update_right_right = Z_update_right_right + Weights(p)*Weights(q)*T;
        %     end
        % end
        % % Z_update_right_right = (Dr(n)/2)^2 * Z_update_right_right;
        % 
        % Z_losses(n,n) = Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right;

        %% Non-Self Term
        for m = n+1:Nep

            % Observation points
            xo1 = first_p(m,1);
            yo1 = first_p(m,2);
            zo1 = first_p(m,3);
            xo2 = second_p(m,1); 
            yo2 = second_p(m,2);
            zo2 = second_p(m,3);
            xo3 = third_p(m,1);
            yo3 = third_p(m,2);
            zo3 = third_p(m,3);
            t_left_obs  = [(xo2-xo1), (yo2-yo1), (zo2-zo1)] / Dl(m);
            t_right_obs = [(xo3-xo2), (yo3-yo2), (zo3-zo2)] / Dr(m);

            %  Testing function
            fm  = zeros(2*order,1); 
            dfm = zeros(2*order,1); 

            % Observation coordinates
            XO  = zeros(2*order,1); 
            YO  = zeros(2*order,1); 
            ZO  = zeros(2*order,1); 
            tlq = zeros(3, 2*order);

            % Parametric integration Point
            lq  = zeros(2*order,1); 

            for q = 1:order

                % Left Segment
                lq(q)        = Dl(m)/2*(Xi(q) + 1);
                fm(q)        = lq(q)/Dl(m);
                dfm(q)       = 1/Dl(m);

                XO(q)        = xo1 + lq(q)/Dl(m)*(xo2-xo1);
                YO(q)        = yo1 + lq(q)/Dl(m)*(yo2-yo1);
                ZO(q)        = zo1 + lq(q)/Dl(m)*(zo2-zo1);
                tlq(:,q) = t_left_obs(:);  % constant tangent for left segment


                % Right Segment
                lq(q+order)  = Dr(m)/2*(Xi(q) + 1);
                fm(q+order)  = (Dr(m)-lq(q+order))/Dr(m);
                dfm(q+order) = -1/Dr(m);

                XO(q+order)  = xo2 + lq(q+order)/Dr(m)*(xo3-xo2);
                YO(q+order)  = yo2 + lq(q+order)/Dr(m)*(yo3-yo2);
                ZO(q+order)  = zo2 + lq(q+order)/Dr(m)*(zo3-zo2);
                tlq(:,q+order) = t_right_obs(:);  % constant tangent for right segment

            end

            % One common segment, source first
            if norm([xs2-xo1 ys2-yo1 zs2-zo1]) == 0
                % Left-Left Segment interactions
                Z_update_left_left = 0;
                for p = 1:order
                    for q = 1:order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_left_left = Z_update_left_left + G*T;
                    end
                end
                Z_update_left_left = Dl(n)/2 * Dl(m)/2 * Z_update_left_left;

                % Right-Right Segment interactions
                Z_update_right_right = 0;
                for p = order+1:2*order
                    for q = order+1:2*order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_right_right = Z_update_right_right + G*T;
                    end
                end
                Z_update_right_right = Dr(n)/2 * Dr(m)/2 * Z_update_right_right;

                % Left-Right Segment interactions
                Z_update_left_right = 0;
                for p = 1:order
                    for q = order+1:2*order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_left_right = Z_update_left_right + G*T;
                    end
                end
                Z_update_left_right = Dl(n)/2 * Dr(m)/2 * Z_update_left_right;

                %Right-Left Segment - Singular
                Z_update_right_left = 0;
                for q = 1:order
                    sqrt1 = sqrt(a^2+lq(q)^2);
                    sqrt2 = sqrt(a^2+(lq(q)-Dr(n))^2);
                    logterm = log(abs(lq(q)+sqrt1)/abs(lq(q)-Dr(n) + sqrt2));
                    S1 = (1-lq(q)/Dr(n))*logterm - 1/Dr(n)*sqrt2 + 1/Dr(n)*sqrt1 - 1i*k0*Dr(n)/2;
                    S2 = -1/Dr(n) * (logterm - 1i*k0*Dr(n));
                    Z_update_right_left = Z_update_right_left + Weights(q)*(fm(q)*S1-1/(k0^2)*dfm(q)*S2);
                end
                Z_update_right_left = Dl(m)/2 * Z_update_right_left;
                Z(m,n) = Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right;

            % One common segment, observation first
            elseif norm([xs1-xo2 ys1-yo2 zs1-zo2]) == 0
                % Left-Left Segment interactions
                Z_update_left_left = 0;
                for p = 1:order
                    for q = 1:order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_left_left = Z_update_left_left + G*T;
                    end
                end
                Z_update_left_left = Dl(n)/2 * Dl(m)/2 * Z_update_left_left;

                % Right-Right Segment interactions
                Z_update_right_right = 0;
                for p = order+1:2*order
                    for q = order+1:2*order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_right_right = Z_update_right_right + G*T;
                    end
                end
                Z_update_right_right = Dr(n)/2 * Dr(m)/2 * Z_update_right_right;

                % Right-Left Segment interactions
                Z_update_right_left = 0;
                for p = order+1:2*order
                    for q = 1:order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_right_left = Z_update_right_left + G*T;
                    end
                end
                Z_update_right_left = Dr(n)/2 * Dl(m)/2 * Z_update_right_left;

                % Left-Right Segment - Singular
                Z_update_left_right = 0;
                for q = order+1:2*order
                    sqrt1 = sqrt(a^2+lq(q)^2);
                    sqrt2 = sqrt(a^2+(lq(q)-Dl(n))^2);
                    logterm = log(abs(lq(q)+sqrt1)/abs(lq(q)-Dl(n) + sqrt2));
                    S1 = lq(q)/Dl(n)*logterm + 1/Dl(n)*sqrt2 -1/Dl(n)*sqrt1 -1i*k0*Dl(n)/2;
                    S2 = 1/Dl(n)*(logterm - 1i*k0*Dl(n));
                    Z_update_left_right = Z_update_left_right + Weights(q)*(fm(q)*S1-1/(k0^2)*dfm(q)*S2);
                end
                Z_update_left_right = Dr(m)/2 * Z_update_left_right;
                Z(m,n) = Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right;

            % No common Element
            else
                % Left-Left Segment interactions
                Z_update_left_left = 0;
                for p = 1:order
                    for q = 1:order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_left_left = Z_update_left_left + G*T;
                    end
                end
                Z_update_left_left = Dl(n)/2 * Dl(m)/2 * Z_update_left_left;
                % Left-Right Segment interactions
                Z_update_left_right = 0;
                for p = 1:order
                    for q = order+1:2*order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_left_right = Z_update_left_right + G*T;
                    end
                end
                Z_update_left_right = Dl(n)/2 * Dr(m)/2 * Z_update_left_right;
                % Right-Left Segment interactions
                Z_update_right_left = 0;
                for p = order+1:2*order
                    for q = 1:order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_right_left = Z_update_right_left + G*T;
                    end
                end
                Z_update_right_left = Dr(n)/2 * Dl(m)/2 * Z_update_right_left;
                % Right-Right Segment interactions
                Z_update_right_right = 0;
                for p = order+1:2*order
                    for q = order+1:2*order
                        T = fn(p)*fm(q)*dot(tlp(:,p),tlq(:,q)) -1/k0^2 * dfn(p)*dfm(q);
                        R = sqrt((XS(p)-XO(q))^2+(YS(p)-YO(q))^2+(ZS(p)-ZO(q))^2 + a^2);
                        G = Weights(p)*Weights(q)*exp(-1i*k0*R)/R;
                        Z_update_right_right = Z_update_right_right + G*T;
                    end
                end
                Z_update_right_right = Dr(n)/2 * Dr(m)/2 * Z_update_right_right;
                % Add all together
                Z(m,n) = Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right;

            end
        end
    end
    Z = sc_const*(Z + Z.' - diag(diag(Z)));
    % Z_losses2 = emc.rho_s*Z_losses;
    Rw = (1/emc.se_copper)/(pi*(2*a-emc.skin_depth)*emc.skin_depth);
    Z_losses = Rw*diag(sqrt(diag(Dl*Dl.')));
    Z = Z + Z_losses; 

end