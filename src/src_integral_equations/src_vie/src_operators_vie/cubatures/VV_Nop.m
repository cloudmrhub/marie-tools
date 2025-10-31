function I = VV_Nop(W, X, Y, Z, Xp, Yp, Zp, k0, r_m, r_n, dx, pwx)

    J = (dx / 2)^6;

    % Compute grids
    Xg = r_m(1) + dx / 2 * X;
    Yg = r_m(2) + dx / 2 * Y;
    Zg = r_m(3) + dx / 2 * Z;

    Xpg = r_n(1) + dx / 2 * Xp;
    Ypg = r_n(2) + dx / 2 * Yp;
    Zpg = r_n(3) + dx / 2 * Zp;

    % Compute GN
    GN = G_N(Xg, Yg, Zg, Xpg, Ypg, Zpg, k0); % Output: [numel(X), 6]

    % Prepare terms for T and B
    Xt = dx / 2 * X;
    Yt = dx / 2 * Y;
    Zt = dx / 2 * Z;

    Xpb = dx / 2 * Xp;
    Ypb = dx / 2 * Yp;
    Zpb = dx / 2 * Zp;

    % Index mapping
    ind_l = [1 2 3 4 1 1 1 2 2 3];
    ind_lp = [1 2 3 4 2 3 4 3 4 4];

    % Precompute T and B for all indices
    T = [ones(size(X)), Xt / dx, Yt / dx, Zt / dx];
    B = [ones(size(Xp)), Xpb / dx, Ypb / dx, Zpb / dx];

    % Preallocate output
    I = zeros(pwx, 6);

    % Vectorize the inner loop
    for llp = 1:pwx
        l = ind_l(llp);
        lp = ind_lp(llp);

        % Compute TB for the current indices
        TB = T(:, l) .* B(:, lp); % Element-wise multiplication

        % Compute the result for all six columns of GN simultaneously
        kernel = (W .* TB) .* GN; % Element-wise multiplication for all columns of GN
        I(llp, :) = sum(kernel, 1) * J; % Sum over rows and scale
    end
end

% function  I = VV_Nop(W,X,Y,Z,Xp,Yp,Zp,k0,r_m,r_n,dx,pwx)
% 
%     J = (dx/2)^6;
% 
%     Xg = r_m(1) + dx/2*X;
%     Yg = r_m(2) + dx/2*Y;
%     Zg = r_m(3) + dx/2*Z;
% 
%     Xpg = r_n(1) + dx/2*Xp;
%     Ypg = r_n(2) + dx/2*Yp;
%     Zpg = r_n(3) + dx/2*Zp;
% 
%     GN = G_N(Xg,Yg,Zg,Xpg,Ypg,Zpg,k0);
% 
%     Xt  = dx/2*X;
%     Yt  = dx/2*Y;
%     Zt  = dx/2*Z;
% 
%     Xpb  = dx/2*Xp;
%     Ypb  = dx/2*Yp;
%     Zpb  = dx/2*Zp;
% 
%     I = zeros(pwx,6);
% 
%     ind_l  = [1 2 3 4 1 1 1 2 2 3];
%     ind_lp = [1 2 3 4 2 3 4 3 4 4];
% 
%     for llp = 1:pwx
% 
%          l = ind_l(llp);
%          lp = ind_lp(llp);
% 
%          switch l
%              case 1
%                  T = 1;
%              case 2
%                  T = Xt/dx;
%              case 3
%                  T = Yt/dx;
%              case 4
%                  T = Zt/dx;
%          end
% 
%          switch lp
%              case 1
%                  B = 1;
%              case 2
%                  B = Xpb/dx;
%              case 3
%                  B = Ypb/dx;
%              case 4
%                  B = Zpb/dx;
%          end
% 
%          TB = T.*B;
% 
%          for i = 1:6
% 
%              kernel = TB .* GN(:,i);
%              I(llp,i) = sum ( W .* kernel ) * J;
% 
%          end
% 
%      end
% 
% end