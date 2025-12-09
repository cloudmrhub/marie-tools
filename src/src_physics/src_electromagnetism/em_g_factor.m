function [g_ax, g_cor, g_sag] = em_g_factor(Phi, b1, cut_ax, cut_cor, cut_sag, Rf, Rp)
% b1:  n1 x n2 x n3 x nc  (B1-)
% Phi: nc x nc  (noise covariance, possibly not PSD)
% Rf, Rp: accelerations along 1st and 2nd in-plane directions

    [n1,n2,n3,nc] = size(b1);

    % Symmetrize Phi and project to PSD
    Phi = (Phi + Phi')/2;
    [V,D] = eig(real(Phi));
    d = diag(D);
    d(d < 0) = 0;
    Phi_psd = V*diag(d)*V';
    Phi_psd = Phi_psd + 1e-12*eye(nc);  % small jitter
    L = chol(Phi_psd,'lower');

    % Axial: x-y plane at z = cut_ax
    if ~isempty(cut_ax) && cut_ax >= 1 && cut_ax <= n3
        b1_ax = squeeze(b1(:,:,cut_ax,:));     % n1 x n2 x nc
        g_ax  = local_g_slice(b1_ax, L, Rf, Rp);
    else
        g_ax = [];
    end

    % Sagittal: y-z plane at x = cut_sag
    if ~isempty(cut_sag) && cut_sag >= 1 && cut_sag <= n1
        b1_sag = squeeze(b1(cut_sag,:,:,:));   % n2 x n3 x nc
        g_sag  = local_g_slice(b1_sag, L, Rf, Rp);
    else
        g_sag = [];
    end

    % Coronal: x-z plane at y = cut_cor
    if ~isempty(cut_cor) && cut_cor >= 1 && cut_cor <= n2
        b1_cor = squeeze(b1(:,cut_cor,:,:));   % n1 x n3 x nc
        g_cor  = local_g_slice(b1_cor, L, Rf, Rp);
    else
        g_cor = [];
    end
end

function g = local_g_slice(b1_slice, L, Rf, Rp)
% b1_slice: nf x np x nc  for a single slice

    [nf,np,nc] = size(b1_slice);

    if mod(nf,Rf) ~= 0 || mod(np,Rp) ~= 0
        error('nf and np must be divisible by Rf and Rp.');
    end

    % --- Whitening ---
    tmp = permute(b1_slice, [3 2 1]);     
    tmp = reshape(tmp, nc, np*nf);        
    tmp = L \ tmp;                        
    tmp = reshape(tmp, nc, np, nf);       
    b1w = permute(tmp, [3 2 1]);          

    Nf_red = nf / Rf;
    Np_red = np / Rp;
    Rtot   = Rf * Rp;

    g = zeros(nf,np);

    for qx = 0:(Nf_red-1)
        full_x = qx + 1 + (0:Rf-1)*Nf_red;

        for qy = 0:(Np_red-1)
            full_y = qy + 1 + (0:Rp-1)*Np_red;

            Rloc = numel(full_x) * numel(full_y);
            S = zeros(nc, Rloc);
            coords = zeros(Rloc,2);
            t = 0;

            for kx = 1:numel(full_x)
                ix = full_x(kx);
                for ky = 1:numel(full_y)
                    iy = full_y(ky);
                    t = t + 1;
                    S(:,t) = squeeze(b1w(ix,iy,:)).';
                    coords(t,:) = [ix,iy];
                end
            end

            % ----- FIX: STABLE INVERSION -----
            A = S' * S;
            eps_reg = 1e-10 * trace(A) / size(A,1);
            Areg = A + eps_reg*eye(size(A));
            Ainv = Areg \ eye(size(A));

            for j = 1:Rloc
                g_val = sqrt(real(Ainv(j,j) * Areg(j,j) / Rtot));
                ix = coords(j,1);
                iy = coords(j,2);
                g(ix,iy) = g_val;
            end
        end
    end
end

