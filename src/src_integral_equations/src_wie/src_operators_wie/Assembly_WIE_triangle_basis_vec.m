function [Z,Z_losses] = Assembly_WIE_triangle_basis_vec(wire,a,order,emc)
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
    Z_diag = zeros(Nep, 1);              % Diagonal terms (Z(n,n))
    Z_lower = cell(Nep, 1);              % Each cell holds [m, Z(m,n)] pairs

    % Gauss quadrature integration
    [Weights,Xi]   = gauss_1d(order);
    Weights        = [Weights;Weights];
    xi_shifted     = (Xi + 1); 
    [pp_ll, qq_ll] = ndgrid(1:order, 1:order);
    [pp_lr, qq_lr] = ndgrid(1:order, order+1:2*order);
    [pp_rl, qq_rl] = ndgrid(order+1:2*order, 1:order);
    [pp_rr, qq_rr] = ndgrid(order+1:2*order, order+1:2*order);

    Weights_ll = Weights(pp_ll) .* Weights(qq_ll);
    Weights_lr = Weights(pp_lr) .* Weights(qq_lr);
    Weights_rl = Weights(pp_rl) .* Weights(qq_rl);
    Weights_rr = Weights(pp_rr) .* Weights(qq_rr);

    % Loop around all points
    parfor n = 1:Nep

        f_p = first_p;
        s_p = second_p;
        t_p = third_p;
        Dl_local = Dl;
        Dr_local = Dr;
        Weights_local_ll = Weights_ll;
        Weights_local_lr = Weights_lr;
        Weights_local_rl = Weights_rl;
        Weights_local_rr = Weights_rr;
        Weights_local = Weights;

        % Source points
        xs1 = f_p(n,1);
        ys1 = f_p(n,2);
        zs1 = f_p(n,3);
        xs2 = s_p(n,1); 
        ys2 = s_p(n,2);
        zs2 = s_p(n,3);
        xs3 = t_p(n,1);
        ys3 = t_p(n,2);
        zs3 = t_p(n,3);

        % Define constant tangent vectors (the segments are straight)
        t_left  = [(xs2-xs1), (ys2-ys1), (zs2-zs1)] / Dl_local(n);
        t_right = [(xs3-xs2), (ys3-ys2), (zs3-zs2)] / Dr_local(n);
        
        % Left segment
        lp_left   = (Dl_local(n)/2) * xi_shifted;                    
        fn_left   = lp_left / Dl_local(n);                          
        dfn_left  = ones(order,1) / Dl_local(n);                    
        XS_left   = xs1 + (lp_left / Dl_local(n)) * (xs2 - xs1);     
        YS_left   = ys1 + (lp_left / Dl_local(n)) * (ys2 - ys1);
        ZS_left   = zs1 + (lp_left / Dl_local(n)) * (zs2 - zs1);
        tlp_left  = repmat(t_left(:), 1, order);               
        
        % Right segment
        lp_right  = (Dr_local(n)/2) * xi_shifted;
        fn_right  = (Dr_local(n) - lp_right) / Dr_local(n);
        dfn_right = -ones(order,1) / Dr_local(n);
        XS_right  = xs2 + (lp_right / Dr_local(n)) * (xs3 - xs2);
        YS_right  = ys2 + (lp_right / Dr_local(n)) * (ys3 - ys2);
        ZS_right  = zs2 + (lp_right / Dr_local(n)) * (zs3 - zs2);
        tlp_right = repmat(t_right(:), 1, order);              
        
        % Combine segments
        lp  = [lp_left; lp_right];               
        fn  = [fn_left; fn_right];
        dfn = [dfn_left; dfn_right];
        XS  = [XS_left; XS_right];
        YS  = [YS_left; YS_right];
        ZS  = [ZS_left; ZS_right];
        tlp = [tlp_left, tlp_right]; 

        %% Self Term
        % Left-Left Segment
        sqrt1   = sqrt(a^2 + (lp(1:order) - Dl_local(n)).^2);
        sqrt2   = sqrt(a^2 + lp(1:order).^2);
        logterm = log(abs(lp(1:order) + sqrt2) ./ abs(lp(1:order) - Dl_local(n) + sqrt1));
        term    = 1i * k0 * Dl_local(n) / 2;
        S1 = (lp(1:order) ./ Dl_local(n)) .* logterm + sqrt1 / Dl_local(n) - sqrt2 / Dl_local(n) - term;
        S2 = (logterm - 2 * term) / Dl_local(n);
        Z_update_left_left = sum(Weights_local(1:order) .* (fn(1:order) .* S1 - (1 / k0^2) * dfn(1:order) .* S2));
        Z_update_left_left = Dl_local(n)/2 * Z_update_left_left;

        % Left-Right Segment               
        dot_tlp = sum(tlp(:, pp_lr) .* tlp(:, qq_lr), 1);      
        dot_tlp = reshape(dot_tlp, order, order);         
        T = (fn(pp_lr) .* fn(qq_lr)) .* dot_tlp - (1/k0^2) * (dfn(pp_lr) .* dfn(qq_lr));         
        dx = XS(pp_lr) - XS(qq_lr);
        dy = YS(pp_lr) - YS(qq_lr);
        dz = ZS(pp_lr) - ZS(qq_lr);
        R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);           
        G = Weights_local_lr .* exp(-1i * k0 * R) ./ R;
        Z_update_left_right = Dl_local(n)/2 * Dr_local(n)/2 * sum(G(:) .* T(:));

        % Right-Left Segment interactions
        dot_tlp = sum(tlp(:, pp_rl) .* tlp(:, qq_rl), 1);
        dot_tlp = reshape(dot_tlp, order, order);
        T = (fn(pp_rl) .* fn(qq_rl)) .* dot_tlp - (1/k0^2) * (dfn(pp_rl) .* dfn(qq_rl));
        dx = XS(pp_rl) - XS(qq_rl);
        dy = YS(pp_rl) - YS(qq_rl);
        dz = ZS(pp_rl) - ZS(qq_rl);
        R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
        G = Weights_local_rl .* exp(-1i * k0 * R) ./ R;
        Z_update_right_left = Dr_local(n)/2 * Dl_local(n)/2 * sum(G(:) .* T(:));

        % Right-Right Segment
        sqrt1   = sqrt(a^2 + (lp(order+1:2*order) - Dr_local(n)).^2);
        sqrt2   = sqrt(a^2 + lp(order+1:2*order).^2);
        logterm = log(abs(lp(order+1:2*order) + sqrt2) ./ abs(lp(order+1:2*order) - Dr_local(n) + sqrt1));
        term    = 1i * k0 * Dr_local(n) / 2;
        S1 = (1 - lp(order+1:2*order) ./ Dr_local(n)) .* logterm - sqrt1 / Dr_local(n) + sqrt2 / Dr_local(n) - term;
        S2 = -(logterm - 2 * term) / Dr_local(n);
        Z_update_right_right = sum(Weights_local(order+1:2*order) .* (fn(order+1:2*order) .* S1 - (1 / k0^2) * dfn(order+1:2*order) .* S2));
        Z_update_right_right = Dr_local(n)/2 * Z_update_right_right;

        % Add all together
        Z_diag(n) = Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right;
        local_Zmn = [];  % Collect [m, val] pairs for Z(m,n)

        %% Non-Self Term
        for m = n+1:Nep

            % Observation points
            xo1 = f_p(m,1);
            yo1 = f_p(m,2);
            zo1 = f_p(m,3);
            xo2 = s_p(m,1); 
            yo2 = s_p(m,2);
            zo2 = s_p(m,3);
            xo3 = t_p(m,1);
            yo3 = t_p(m,2);
            zo3 = t_p(m,3);
            t_left_obs  = [(xo2-xo1), (yo2-yo1), (zo2-zo1)] / Dl_local(m);
            t_right_obs = [(xo3-xo2), (yo3-yo2), (zo3-zo2)] / Dr_local(m);

            % Left segment 
            lq_left   = (Dl_local(m)/2) * xi_shifted;                          
            fm_left   = lq_left / Dl_local(m);                                 
            dfm_left  = ones(order,1) / Dl_local(m);                           
            XO_left   = xo1 + (lq_left / Dl_local(m)) * (xo2 - xo1);
            YO_left   = yo1 + (lq_left / Dl_local(m)) * (yo2 - yo1);
            ZO_left   = zo1 + (lq_left / Dl_local(m)) * (zo2 - zo1);
            tlq_left  = repmat(t_left_obs(:), 1, order);                 
            
            % Right segment
            lq_right  = (Dr_local(m)/2) * xi_shifted;
            fm_right  = (Dr_local(m) - lq_right) / Dr_local(m);
            dfm_right = -ones(order,1) / Dr_local(m);
            XO_right  = xo2 + (lq_right / Dr_local(m)) * (xo3 - xo2);
            YO_right  = yo2 + (lq_right / Dr_local(m)) * (yo3 - yo2);
            ZO_right  = zo2 + (lq_right / Dr_local(m)) * (zo3 - zo2);
            tlq_right = repmat(t_right_obs(:), 1, order);              
            
            % Combine segments
            lq  = [lq_left; lq_right];             
            fm  = [fm_left; fm_right];             
            dfm = [dfm_left; dfm_right];
            XO  = [XO_left; XO_right];             
            YO  = [YO_left; YO_right];
            ZO  = [ZO_left; ZO_right];
            tlq = [tlq_left, tlq_right];           

            % One common segment, source first
            if norm([xs2-xo1 ys2-yo1 zs2-zo1]) == 0
                % Left-Left Segment interactions
                dot_tlpq = sum(tlp(:, pp_ll) .* tlq(:, qq_ll), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_ll) .* fm(qq_ll)) .* dot_tlpq - (1/k0^2) * (dfn(pp_ll) .* dfm(qq_ll));
                dx = XS(pp_ll) - XO(qq_ll);
                dy = YS(pp_ll) - YO(qq_ll);
                dz = ZS(pp_ll) - ZO(qq_ll);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_ll .* exp(-1i * k0 * R) ./ R;
                Z_update_left_left = Dl_local(n)/2 * Dl_local(m)/2 * sum(G(:) .* T(:));

                % Right-Right Segment interactions
                dot_tlpq = sum(tlp(:, pp_rr) .* tlq(:, qq_rr), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_rr) .* fm(qq_rr)) .* dot_tlpq - (1/k0^2) * (dfn(pp_rr) .* dfm(qq_rr));
                dx = XS(pp_rr) - XO(qq_rr);
                dy = YS(pp_rr) - YO(qq_rr);
                dz = ZS(pp_rr) - ZO(qq_rr);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_rr .* exp(-1i * k0 * R) ./ R;
                Z_update_right_right = Dr_local(n)/2 * Dr_local(m)/2 * sum(G(:) .* T(:));

                % Left-Right Segment interactions
                dot_tlpq = sum(tlp(:, pp_lr) .* tlq(:, qq_lr), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_lr) .* fm(qq_lr)) .* dot_tlpq - (1/k0^2) * (dfn(pp_lr) .* dfm(qq_lr));
                dx = XS(pp_lr) - XO(qq_lr);
                dy = YS(pp_lr) - YO(qq_lr);
                dz = ZS(pp_lr) - ZO(qq_lr);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_lr .* exp(-1i * k0 * R) ./ R;
                Z_update_left_right = Dl_local(n)/2 * Dr_local(m)/2 * sum(G(:) .* T(:));

                %Right-Left Segment - Singular
                sqrt1   = sqrt(a^2 + lq(1:order).^2);
                sqrt2   = sqrt(a^2 + (lq(1:order) - Dr_local(n)).^2);
                logterm = log(abs(lq(1:order) + sqrt1) ./ abs(lq(1:order) - Dr_local(n) + sqrt2));
                term    = 1i * k0 * Dr_local(n) / 2;
                S1 = (1 - lq(1:order) / Dr_local(n)) .* logterm - sqrt2 / Dr_local(n) + sqrt1 / Dr_local(n) - term;
                S2 = -(logterm - term * 2) / Dr_local(n);
                Z_update_right_left = sum(Weights_local(1:order) .* (fm(1:order) .* S1 - (1 / k0^2) * dfm(1:order) .* S2));
                Z_update_right_left = Dl_local(m)/2 * Z_update_right_left;

                % Add all together
                local_Zmn = [local_Zmn; m, Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right];

            % One common segment, observation first
            elseif norm([xs1-xo2 ys1-yo2 zs1-zo2]) == 0
                % Left-Left Segment interactions
                dot_tlpq = sum(tlp(:, pp_ll) .* tlq(:, qq_ll), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_ll) .* fm(qq_ll)) .* dot_tlpq - (1 / k0^2) * (dfn(pp_ll) .* dfm(qq_ll));
                dx = XS(pp_ll) - XO(qq_ll);
                dy = YS(pp_ll) - YO(qq_ll);
                dz = ZS(pp_ll) - ZO(qq_ll);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_ll .* exp(-1i * k0 * R) ./ R;
                Z_update_left_left = Dl_local(n)/2 * Dl_local(m)/2 * sum(G(:) .* T(:));

                % Right-Right Segment interactions
                dot_tlpq = sum(tlp(:, pp_rr) .* tlq(:, qq_rr), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_rr) .* fm(qq_rr)) .* dot_tlpq - (1 / k0^2) * (dfn(pp_rr) .* dfm(qq_rr));
                dx = XS(pp_rr) - XO(qq_rr);
                dy = YS(pp_rr) - YO(qq_rr);
                dz = ZS(pp_rr) - ZO(qq_rr);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_rr .* exp(-1i * k0 * R) ./ R;
                Z_update_right_right = Dr_local(n)/2 * Dr_local(m)/2 * sum(G(:) .* T(:));

                % Right-Left Segment interactions
                dot_tlpq = sum(tlp(:, pp_rl) .* tlq(:, qq_rl), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_rl) .* fm(qq_rl)) .* dot_tlpq - (1 / k0^2) * (dfn(pp_rl) .* dfm(qq_rl));
                dx = XS(pp_rl) - XO(qq_rl);
                dy = YS(pp_rl) - YO(qq_rl);
                dz = ZS(pp_rl) - ZO(qq_rl);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_rl .* exp(-1i * k0 * R) ./ R;
                Z_update_right_left = Dr_local(n)/2 * Dl_local(m)/2 * sum(G(:) .* T(:));

                % Left-Right Segment - Singular
                sqrt1   = sqrt(a^2 + lq(order+1:2*order).^2);
                sqrt2   = sqrt(a^2 + (lq(order+1:2*order) - Dl_local(n)).^2);
                logterm = log(abs(lq(order+1:2*order) + sqrt1) ./ abs(lq(order+1:2*order) - Dl_local(n) + sqrt2));
                term    = 1i * k0 * Dl_local(n) / 2;
                S1 = (lq(order+1:2*order) / Dl_local(n)) .* logterm + sqrt2 / Dl_local(n) - sqrt1 / Dl_local(n) - term;
                S2 = (logterm - 2 * term) / Dl_local(n);
                Z_update_left_right = sum(Weights_local(order+1:2*order) .* (fm(order+1:2*order) .* S1 - (1 / k0^2) * dfm(order+1:2*order) .* S2));
                Z_update_left_right = Dr_local(m)/2 * Z_update_left_right;

                % Add all together
                local_Zmn = [local_Zmn; m, Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right];

            % No common Element
            else
                % Left-Left Segment interactions
                dot_tlpq = sum(tlp(:, pp_ll) .* tlq(:, qq_ll), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_ll) .* fm(qq_ll)) .* dot_tlpq - (1 / k0^2) * (dfn(pp_ll) .* dfm(qq_ll));
                dx = XS(pp_ll) - XO(qq_ll);
                dy = YS(pp_ll) - YO(qq_ll);
                dz = ZS(pp_ll) - ZO(qq_ll);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_ll .* exp(-1i * k0 * R) ./ R;
                Z_update_left_left = Dl_local(n)/2 * Dl_local(m)/2 * sum(G(:) .* T(:));

                % Left-Right Segment interactions
                dot_tlpq = sum(tlp(:, pp_lr) .* tlq(:, qq_lr), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_lr) .* fm(qq_lr)) .* dot_tlpq - (1 / k0^2) * (dfn(pp_lr) .* dfm(qq_lr));
                dx = XS(pp_lr) - XO(qq_lr);
                dy = YS(pp_lr) - YO(qq_lr);
                dz = ZS(pp_lr) - ZO(qq_lr);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_lr .* exp(-1i * k0 * R) ./ R;
                Z_update_left_right = Dl_local(n)/2 * Dr_local(m)/2 * sum(G(:) .* T(:));

                % Right-Left Segment interactions
                dot_tlpq = sum(tlp(:, pp_rl) .* tlq(:, qq_rl), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_rl) .* fm(qq_rl)) .* dot_tlpq - (1 / k0^2) * (dfn(pp_rl) .* dfm(qq_rl));
                dx = XS(pp_rl) - XO(qq_rl);
                dy = YS(pp_rl) - YO(qq_rl);
                dz = ZS(pp_rl) - ZO(qq_rl);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_rl .* exp(-1i * k0 * R) ./ R;
                Z_update_right_left = Dr_local(n)/2 * Dl_local(m)/2 * sum(G(:) .* T(:));

                % Right-Right Segment interactions
                dot_tlpq = sum(tlp(:, pp_rr) .* tlq(:, qq_rr), 1);
                dot_tlpq = reshape(dot_tlpq, order, order);
                T = (fn(pp_rr) .* fm(qq_rr)) .* dot_tlpq - (1 / k0^2) * (dfn(pp_rr) .* dfm(qq_rr));
                dx = XS(pp_rr) - XO(qq_rr);
                dy = YS(pp_rr) - YO(qq_rr);
                dz = ZS(pp_rr) - ZO(qq_rr);
                R = sqrt(dx.^2 + dy.^2 + dz.^2 + a^2);
                G = Weights_local_rr .* exp(-1i * k0 * R) ./ R;
                Z_update_right_right = Dr_local(n)/2 * Dr_local(m)/2 * sum(G(:) .* T(:));

                % Add all together
                local_Zmn = [local_Zmn; m, Z_update_left_left + Z_update_left_right + Z_update_right_left + Z_update_right_right];

            end
            Z_lower{n} = local_Zmn;  % Store result for column n
        end
    end

    Z = zeros(Nep);
    for n = 1:Nep
        Z(n,n) = Z_diag(n);
        if ~isempty(Z_lower{n})
            rows = Z_lower{n}(:,1);
            vals = Z_lower{n}(:,2);
            Z(rows, n) = vals;
        end
    end
    
    Z = sc_const*(Z + Z.' - diag(diag(Z)));
    Rw = (1/emc.se_copper)/(pi*(2*a-emc.skin_depth)*emc.skin_depth);
    % Z_losses = Rw*diag(sqrt(diag(Dl*Dl.')));
    Z_losses = spdiags(Rw*abs(Dl), 0, n, n);
    Z = Z + Z_losses; 

end