function [P,ndeim,xds,yds,zds] = rSVD_deim(U,ql,Nscat,idxS,r)

    % get sizes of the input orthonormal matrix
    [n,m] = size(U);

    % phi is the vector with the indexes
    phi = zeros(m, 1);

    % first index: maximum value
    [~, phi(1)] = max(abs(U(:, 1)));

    % loop on the columns of U
    for j = 2:m
        idx = phi(1:j-1);
        uL = U(:,j); % select column vector
        c = uL(idx); % select entries of current column vector
        c = U(idx,1:j-1)\c; % compute coeff
        resid = uL - U(:,1:j-1)*c; % orthogonalize
        [~, phi(j)] = max(abs(resid)); % find maximum entry in column

    end
    % prepare elements for return
    phi = sort(phi);
    e = speye(n);
    P = sparse(e(:, phi));

    xs = r(:,:,:,1);
    ys = r(:,:,:,2);
    zs = r(:,:,:,3);
    xs = xs(idxS);
    ys = ys(idxS);
    zs = zs(idxS);
    clear r % no longer needed
    if ql == 3
        P3D = [P(1:Nscat,:), P(Nscat+1:2*Nscat,:), P(2*Nscat+1:3*Nscat,:)];
    elseif ql == 12
        P3D = [P(1:Nscat,:), P(Nscat+1:2*Nscat,:), P(2*Nscat+1:3*Nscat,:), P(3*Nscat+1:4*Nscat,:), P(4*Nscat+1:5*Nscat,:), P(5*Nscat+1:6*Nscat,:), P(6*Nscat+1:7*Nscat,:), P(7*Nscat+1:8*Nscat,:), P(8*Nscat+1:9*Nscat,:), P(9*Nscat+1:10*Nscat,:), P(10*Nscat+1:11*Nscat,:), P(11*Nscat+1:12*Nscat,:)];
    end
    V3D = sum(P3D,2);
    clear P3D  % no longer needed
    idxD = find(V3D);
    idxD = sort(idxD);
    xds = xs(idxD);
    yds = ys(idxD);
    zds = zs(idxD);
    clear xs ys zs V3D  % no longer needed
    if ql == 3
        idxD = [idxD; Nscat+idxD; 2*Nscat+idxD];
        e = speye(3*Nscat);
    elseif ql == 12
        idxD = [idxD; Nscat+idxD; 2*Nscat+idxD; 3*Nscat+idxD; 4*Nscat+idxD; 5*Nscat+idxD; 6*Nscat+idxD; 7*Nscat+idxD; 8*Nscat+idxD; 9*Nscat+idxD; 10*Nscat+idxD; 11*Nscat+idxD];
        e = speye(12*Nscat);
    end
    P = sparse(e(:,idxD));
    ndeim = size(P,2);

end